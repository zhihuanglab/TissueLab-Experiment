#!/usr/bin/env python3
"""Paired held-out comparisons of the SEA-AD bulk biomarker and WSI baselines."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


OUTCOME = "slope_zmem0"
COVARIATES = ["max_age_vis", "sex", "braak_numeric", "cerad_ordinal"]
BULK_COLUMN = "round_0002__ca1_reactive_astrocyte_lineage_fraction"

BASELINES = [
    ("UNI-ABMIL", "uni_gated_abmil/test_predictions.csv", "model_prediction"),
    ("MUSK-ABMIL", "musk_gated_abmil/test_predictions.csv", "model_prediction"),
    ("TITAN Linear Probe", "titan_linear_probe/test_predictions.csv", "model_prediction"),
    ("TITAN Fine-tune", "titan_finetune/test_predictions.csv", "titan_prediction"),
]


def _sex_numeric(series: pd.Series) -> np.ndarray:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map({"female": 0.0, "f": 0.0, "male": 1.0, "m": 1.0})
        .to_numpy(dtype=float)
    )


def _covariate_matrix(frame: pd.DataFrame) -> np.ndarray:
    return np.column_stack(
        [
            np.ones(len(frame)),
            frame["max_age_vis"].to_numpy(dtype=float),
            _sex_numeric(frame["sex"]),
            frame["braak_numeric"].to_numpy(dtype=float),
            frame["cerad_ordinal"].to_numpy(dtype=float),
        ]
    )


def _residualize(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    return values - covariates @ np.linalg.lstsq(covariates, values, rcond=None)[0]


def _partial_r(frame: pd.DataFrame, column: str, direction: int) -> float:
    x = frame[column].to_numpy(dtype=float)
    y = frame[OUTCOME].to_numpy(dtype=float)
    z = _covariate_matrix(frame)
    x_residual = _residualize(x, z)
    y_residual = _residualize(y, z)
    if x_residual.std() <= 1e-15 or y_residual.std() <= 1e-15:
        return float("nan")
    return float(direction * np.corrcoef(x_residual, y_residual)[0, 1])


def _bootstrap(
    frame: pd.DataFrame,
    methods: dict[str, tuple[str, int]],
    *,
    repetitions: int,
    seed: int,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    method_values = {name: [] for name in methods}
    differences = {name: [] for name in methods if name != "TissueLab bulk CA1 reactivity"}
    for _ in range(repetitions):
        index = rng.integers(0, len(frame), len(frame))
        sample = frame.iloc[index].reset_index(drop=True)
        values = {
            name: _partial_r(sample, column, direction)
            for name, (column, direction) in methods.items()
        }
        for name, value in values.items():
            if np.isfinite(value):
                method_values[name].append(value)
        tissue_value = values["TissueLab bulk CA1 reactivity"]
        for name in differences:
            if np.isfinite(tissue_value) and np.isfinite(values[name]):
                differences[name].append(tissue_value - values[name])
    return (
        {name: np.asarray(values, dtype=float) for name, values in method_values.items()},
        {name: np.asarray(values, dtype=float) for name, values in differences.items()},
    )


def _paired_permutation_p(
    frame: pd.DataFrame,
    tissue_column: str,
    baseline_column: str,
    *,
    repetitions: int,
    seed: int,
) -> tuple[float, float]:
    covariates = _covariate_matrix(frame)
    outcome = _residualize(frame[OUTCOME].to_numpy(float), covariates)
    tissue = -_residualize(frame[tissue_column].to_numpy(float), covariates)
    baseline = _residualize(frame[baseline_column].to_numpy(float), covariates)

    def standardize(values: np.ndarray) -> np.ndarray:
        return (values - values.mean()) / values.std(ddof=1)

    outcome = standardize(outcome)
    tissue = standardize(tissue)
    baseline = standardize(baseline)
    observed = float(np.corrcoef(tissue, outcome)[0, 1] - np.corrcoef(baseline, outcome)[0, 1])

    rng = np.random.default_rng(seed)
    exceed = 0
    completed = 0
    chunk_size = 5000
    while completed < repetitions:
        count = min(chunk_size, repetitions - completed)
        swap = rng.integers(0, 2, size=(count, len(frame)), dtype=np.int8).astype(bool)
        first = np.where(swap, baseline, tissue)
        second = np.where(swap, tissue, baseline)
        first = first - first.mean(axis=1, keepdims=True)
        second = second - second.mean(axis=1, keepdims=True)
        first /= first.std(axis=1, ddof=1, keepdims=True)
        second /= second.std(axis=1, ddof=1, keepdims=True)
        permutation_delta = (
            np.sum(first * outcome, axis=1) / (len(frame) - 1)
            - np.sum(second * outcome, axis=1) / (len(frame) - 1)
        )
        exceed += int(np.sum(np.abs(permutation_delta) >= abs(observed)))
        completed += count
    return observed, float((exceed + 1) / (repetitions + 1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--baseline-dir",
        type=Path,
        default=Path("docs/ad_baseline_standardized_20260826"),
    )
    parser.add_argument(
        "--feature-matrix",
        type=Path,
        default=Path(
            "docs/ad_april_reconstructed_statistics_20260825/inputs/"
            "run_20260420_162449_64fb7a_test.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/ad_baseline_standardized_20260826/comparisons"),
    )
    parser.add_argument("--bootstrap-repetitions", type=int, default=10_000)
    parser.add_argument("--permutation-repetitions", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20260826)
    args = parser.parse_args()

    feature_frame = pd.read_csv(args.feature_matrix)
    frame = feature_frame[
        ["donor_id", OUTCOME, *COVARIATES, BULK_COLUMN]
    ].copy()
    methods: dict[str, tuple[str, int]] = {
        "TissueLab bulk CA1 reactivity": (BULK_COLUMN, -1)
    }
    for name, relative_path, prediction_column in BASELINES:
        predictions = pd.read_csv(args.baseline_dir / relative_path)
        renamed = f"prediction__{name.replace(' ', '_')}"
        frame = frame.merge(
            predictions[["donor_id", prediction_column]].rename(
                columns={prediction_column: renamed}
            ),
            on="donor_id",
            validate="one_to_one",
        )
        methods[name] = (renamed, 1)

    if len(frame) != 34 or frame["donor_id"].duplicated().any():
        raise ValueError("Expected 34 unique held-out donors")

    observed = {
        name: _partial_r(frame, column, direction)
        for name, (column, direction) in methods.items()
    }
    method_bootstrap, difference_bootstrap = _bootstrap(
        frame,
        methods,
        repetitions=args.bootstrap_repetitions,
        seed=args.seed,
    )

    method_rows = []
    for name in methods:
        low, high = np.quantile(method_bootstrap[name], [0.025, 0.975])
        method_rows.append(
            {
                "method": name,
                "n": len(frame),
                "direction_aligned_partial_r": observed[name],
                "bootstrap_ci_low": float(low),
                "bootstrap_ci_high": float(high),
                "bootstrap_repetitions": args.bootstrap_repetitions,
            }
        )

    comparison_rows = []
    raw_p = []
    for offset, (name, _, _) in enumerate(BASELINES):
        values = difference_bootstrap[name]
        low, high = np.quantile(values, [0.025, 0.975])
        permutation_delta, permutation_p = _paired_permutation_p(
            frame,
            BULK_COLUMN,
            methods[name][0],
            repetitions=args.permutation_repetitions,
            seed=args.seed + 10_000 + offset,
        )
        delta = observed["TissueLab bulk CA1 reactivity"] - observed[name]
        if not np.isclose(delta, permutation_delta, atol=1e-12):
            raise RuntimeError(f"Delta mismatch for {name}: {delta} vs {permutation_delta}")
        raw_p.append(permutation_p)
        comparison_rows.append(
            {
                "comparison": f"TissueLab bulk CA1 reactivity minus {name}",
                "n": len(frame),
                "tissuelab_partial_r": observed["TissueLab bulk CA1 reactivity"],
                "baseline_partial_r": observed[name],
                "delta_partial_r": delta,
                "paired_bootstrap_ci_low": float(low),
                "paired_bootstrap_ci_high": float(high),
                "paired_permutation_p_two_sided": permutation_p,
                "bootstrap_repetitions": args.bootstrap_repetitions,
                "permutation_repetitions": args.permutation_repetitions,
            }
        )

    adjusted = multipletests(raw_p, method="holm")[1]
    for row, value in zip(comparison_rows, adjusted):
        row["holm_p_across_four_baselines"] = float(value)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    methods_table = pd.DataFrame(method_rows)
    comparisons_table = pd.DataFrame(comparison_rows)
    methods_table.to_csv(args.output_dir / "heldout_method_partial_r.csv", index=False)
    comparisons_table.to_csv(
        args.output_dir / "tissuelab_vs_wsi_baselines_paired.csv", index=False
    )
    summary = {
        "n": len(frame),
        "outcome": OUTCOME,
        "covariates": COVARIATES,
        "primary_tissuelab_feature": "bulk CA1 reactive-astrocyte fraction",
        "direction_alignment": "TissueLab multiplied by -1; baseline predictions unchanged",
        "bootstrap_repetitions": args.bootstrap_repetitions,
        "permutation_repetitions": args.permutation_repetitions,
        "seed": args.seed,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(methods_table.to_string(index=False))
    print()
    print(comparisons_table.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Shared-donor bootstrap for the fitted SEA-AD panel association."""

from __future__ import annotations

import argparse
import json
import re
from statistics import NormalDist
from pathlib import Path

import numpy as np
import pandas as pd


COVARIATES = ["max_age_vis", "sex_binary", "braak_numeric", "cerad_ordinal"]
ROUND_RE = re.compile(r"^round_(\d{4})__")
NORMAL = NormalDist()


def fitted_multiple_r(
    outcome: np.ndarray,
    covariates: np.ndarray,
    features: np.ndarray,
) -> float:
    design = np.column_stack([np.ones(len(outcome)), covariates])
    residual_maker = np.eye(len(outcome)) - design @ np.linalg.pinv(design)
    outcome_residual = residual_maker @ outcome
    feature_residual = residual_maker @ features
    feature_residual = feature_residual[:, feature_residual.std(axis=0) > 1e-12]
    if feature_residual.shape[1] == 0 or np.linalg.norm(outcome_residual) <= 1e-12:
        return np.nan
    prediction = feature_residual @ np.linalg.lstsq(
        feature_residual, outcome_residual, rcond=None
    )[0]
    if np.std(prediction) <= 1e-12:
        return np.nan
    return float(np.corrcoef(prediction, outcome_residual)[0, 1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-samples", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--round", type=int, default=10)
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("run_*_test.csv"))
    if len(paths) != 5:
        raise ValueError(f"Expected five Test files, found {len(paths)}")

    frames = [
        pd.read_csv(path).dropna(subset=["slope_zmem0"]).reset_index(drop=True)
        for path in paths
    ]
    first = frames[0]
    donor_order = first["donor_id"].astype(str).tolist()
    for path, frame in zip(paths[1:], frames[1:]):
        if frame["donor_id"].astype(str).tolist() != donor_order:
            raise ValueError(f"Donor order differs in {path}")
        if not np.allclose(frame["slope_zmem0"], first["slope_zmem0"]):
            raise ValueError(f"Outcome differs in {path}")
        if not np.allclose(frame[COVARIATES], first[COVARIATES]):
            raise ValueError(f"Covariates differ in {path}")

    outcome = first["slope_zmem0"].to_numpy(float)
    covariates = first[COVARIATES].to_numpy(float)
    features_by_run = []
    for path, frame in zip(paths, frames):
        columns = [
            column
            for column in frame.columns
            if (match := ROUND_RE.match(column))
            and int(match.group(1)) <= args.round
        ]
        if not columns:
            raise ValueError(f"No round-{args.round} features in {path}")
        if frame[columns].isna().any().any():
            raise ValueError(f"Missing round-{args.round} feature values in {path}")
        features_by_run.append(frame[columns].to_numpy(float))

    observed_by_run = np.array(
        [fitted_multiple_r(outcome, covariates, features) for features in features_by_run]
    )
    observed_mean = float(observed_by_run.mean())

    rng = np.random.default_rng(args.seed)
    bootstrap_values = np.empty(args.bootstrap_samples, dtype=float)
    for iteration in range(args.bootstrap_samples):
        indices = rng.integers(0, len(outcome), size=len(outcome))
        run_values = [
            fitted_multiple_r(
                outcome[indices], covariates[indices], features[indices]
            )
            for features in features_by_run
        ]
        bootstrap_values[iteration] = np.nanmean(run_values)

    valid = bootstrap_values[np.isfinite(bootstrap_values)]
    if len(valid) != args.bootstrap_samples:
        raise RuntimeError(
            f"Only {len(valid)} of {args.bootstrap_samples} bootstrap samples were valid"
        )

    percentile_ci = np.percentile(valid, [2.5, 97.5])
    basic_ci = np.array(
        [2 * observed_mean - percentile_ci[1], 2 * observed_mean - percentile_ci[0]]
    )

    jackknife_values = np.empty(len(outcome), dtype=float)
    all_indices = np.arange(len(outcome))
    for omitted in range(len(outcome)):
        indices = all_indices[all_indices != omitted]
        jackknife_values[omitted] = np.mean(
            [
                fitted_multiple_r(
                    outcome[indices], covariates[indices], features[indices]
                )
                for features in features_by_run
            ]
        )
    jackknife_center = jackknife_values.mean()
    jackknife_delta = jackknife_center - jackknife_values
    acceleration = float(
        np.sum(jackknife_delta**3)
        / (6 * np.sum(jackknife_delta**2) ** 1.5)
    )
    less_fraction = float(np.mean(valid < observed_mean))
    less_fraction = min(max(less_fraction, 0.5 / len(valid)), 1 - 0.5 / len(valid))
    bias_correction = NORMAL.inv_cdf(less_fraction)
    adjusted_probabilities = []
    for alpha in [0.025, 0.975]:
        z_alpha = NORMAL.inv_cdf(alpha)
        adjusted_z = bias_correction + (
            (bias_correction + z_alpha)
            / (1 - acceleration * (bias_correction + z_alpha))
        )
        adjusted_probabilities.append(NORMAL.cdf(adjusted_z))
    bca_ci = np.quantile(valid, adjusted_probabilities)

    payload = {
        "statistic": "mean within-Test-cohort fitted multiple partial R across five runs",
        "interpretation": "fitted association, not frozen out-of-sample panel prediction",
        "cohort": "Test",
        "n_donors": len(outcome),
        "round": args.round,
        "bootstrap": "shared-donor nonparametric bootstrap with complete refitting",
        "bootstrap_samples": args.bootstrap_samples,
        "seed": args.seed,
        "observed_by_run": observed_by_run.tolist(),
        "observed_mean": observed_mean,
        "observed_run_sd_ddof1": float(observed_by_run.std(ddof=1)),
        "observed_run_min": float(observed_by_run.min()),
        "observed_run_max": float(observed_by_run.max()),
        "bootstrap_mean": float(valid.mean()),
        "bootstrap_median": float(np.median(valid)),
        "bootstrap_standard_error": float(valid.std(ddof=1)),
        "bootstrap_ci95_percentile": percentile_ci.tolist(),
        "bootstrap_ci95_basic": basic_ci.tolist(),
        "bootstrap_bca_bias_correction": bias_correction,
        "bootstrap_bca_acceleration": acceleration,
        "bootstrap_bca_adjusted_quantiles": adjusted_probabilities,
        "bootstrap_ci95_bca": bca_ci.tolist(),
        "jackknife_min": float(jackknife_values.min()),
        "jackknife_max": float(jackknife_values.max()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

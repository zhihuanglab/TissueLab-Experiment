#!/usr/bin/env python3
"""Confirmatory shared-donor Freedman-Lane test for Test round-10 trajectory."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


COVARIATES = ["max_age_vis", "sex_binary", "braak_numeric", "cerad_ordinal"]
ROUND_RE = re.compile(r"^round_(\d{4})__")


def projection_basis(matrix: np.ndarray, tolerance: float = 1e-12) -> np.ndarray:
    u, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    if len(singular) == 0:
        return np.empty((len(matrix), 0), dtype=float)
    keep = singular > tolerance * singular[0]
    return u[:, keep]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--batch-size", type=int, default=2000)
    parser.add_argument("--cohort", choices=["train", "test"], default="test")
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob(f"run_*_{args.cohort}.csv"))
    if len(paths) != 5:
        raise ValueError(f"Expected five Test files, found {len(paths)}")
    frames = [pd.read_csv(path).dropna(subset=["slope_zmem0"]).reset_index(drop=True) for path in paths]
    donor_order = frames[0]["donor_id"].astype(str).tolist()
    if any(frame["donor_id"].astype(str).tolist() != donor_order for frame in frames[1:]):
        raise ValueError("Test donor order differs across runs")

    first = frames[0]
    design = np.column_stack([np.ones(len(first)), first[COVARIATES].to_numpy(float)])
    outcome = first["slope_zmem0"].to_numpy(float)
    residual_maker = np.eye(len(first)) - design @ np.linalg.pinv(design)
    outcome_residual = residual_maker @ outcome

    bases_by_run = []
    observed_by_run_curve = []
    for frame in frames:
        if not np.allclose(frame["slope_zmem0"].to_numpy(float), outcome):
            raise ValueError("Test outcome differs across runs")
        if not np.allclose(frame[COVARIATES].to_numpy(float), first[COVARIATES].to_numpy(float)):
            raise ValueError("Test covariates differ across runs")
        run_bases = []
        run_observed = []
        for round_id in range(1, 11):
            columns = [
                column for column in frame.columns
                if (match := ROUND_RE.match(column)) and int(match.group(1)) <= round_id
            ]
            if frame[columns].isna().any().any():
                raise ValueError(f"Missing round-{round_id} feature values in {paths[len(bases_by_run)]}")
            matrix = residual_maker @ frame[columns].to_numpy(float)
            matrix = matrix[:, matrix.std(axis=0) > 1e-12]
            basis = projection_basis(matrix)
            run_bases.append(basis)
            run_observed.append(float(np.linalg.norm(basis.T @ outcome_residual) / np.linalg.norm(outcome_residual)))
        bases_by_run.append(run_bases)
        observed_by_run_curve.append(run_observed)

    observed_curve = np.mean(observed_by_run_curve, axis=0)
    observed = float(observed_curve[-1])
    rng = np.random.default_rng(args.seed)
    exceedances = np.zeros(10, dtype=np.int64)
    null_values = np.empty((args.permutations, 10), dtype=np.float32)
    for start in range(0, args.permutations, args.batch_size):
        end = min(args.permutations, start + args.batch_size)
        permuted = np.column_stack([
            outcome_residual[rng.permutation(len(outcome_residual))]
            for _ in range(end - start)
        ])
        residualized = residual_maker @ permuted
        denominator = np.linalg.norm(residualized, axis=0)
        for round_index in range(10):
            run_statistics = np.vstack([
                np.linalg.norm(run_bases[round_index].T @ residualized, axis=0) / denominator
                for run_bases in bases_by_run
            ])
            mean_statistics = run_statistics.mean(axis=0)
            null_values[start:end, round_index] = mean_statistics
            exceedances[round_index] += int(np.sum(mean_statistics >= observed_curve[round_index]))

    payload = {
        "cohort": args.cohort.title(),
        "round": 10,
        "statistic": "mean within-cohort fitted partial correlation across five runs",
        "permutation": "shared-donor Freedman-Lane residual permutation",
        "permutations": args.permutations,
        "seed": args.seed,
        "n_donors": len(outcome),
        "observed_by_run": [row[-1] for row in observed_by_run_curve],
        "observed_mean": observed,
        "null_median": float(np.percentile(null_values[:, -1], 50)),
        "null_p95": float(np.percentile(null_values[:, -1], 95)),
        "null_p99": float(np.percentile(null_values[:, -1], 99)),
        "exceedances": int(exceedances[-1]),
        "permutation_p_upper": float((1 + exceedances[-1]) / (1 + args.permutations)),
        "trajectory": {
            "observed": observed_curve.tolist(),
            "null_median": np.percentile(null_values, 50, axis=0).tolist(),
            "null_p95": np.percentile(null_values, 95, axis=0).tolist(),
            "permutation_p_upper": ((1 + exceedances) / (1 + args.permutations)).tolist(),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

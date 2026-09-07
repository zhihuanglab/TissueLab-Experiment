#!/usr/bin/env python3
"""Held-out statistics for the original April SEA-AD measurement stack."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


RUNS = {
    "run1": "run_20260420_162449_64fb7a",
    "run2": "run_20260420_215745_90b105",
    "run3": "run_20260421_001115_585b70",
    "run4": "run_20260421_014751_94b13b",
    "run5": "run_20260421_052101_efbee1",
}

BULK_SPEC = (
    "run1",
    "round_0002__ca1_reactive_astrocyte_lineage_fraction",
    -1,
    "CA1 reactive-astrocyte fraction",
)

SPATIAL_SPECS = [
    ("run1", "round_0003__ca1_pyramidal_reactive_astro_niche_fraction_80px", -1, ">=1 reactive astrocyte near CA1 pyramidal neuron, approximately 40 um"),
    ("run1", "round_0004__ca1_pyramidal_reactive_astro_coverage_60px", -1, ">=1 reactive astrocyte near CA1 pyramidal neuron, approximately 30 um"),
    ("run2", "round_0003__ca1_pyramidal_reactive_niche__reactive_over_astroglial_neighbored_pyramidal_fraction_r70px", -1, "Reactive-neighbored / astroglial-neighbored pyramidal ratio, approximately 35 um"),
    ("run3", "round_0002__peripyramidal_reactivity_r35um", -1, "Peripyramidal astrocyte-lineage reactivity, 35 um"),
    ("run3", "round_0003__ca1_pyramidal_cuffed_fraction_r35um_ge1", -1, ">=1 reactive astrocyte within 35 um"),
    ("run3", "round_0004__ca1_pyramidal_cuffed_fraction_r35um_ge2", -1, ">=2 reactive astrocytes within 35 um"),
    ("run3", "round_0005__ca1_pyramidal_supermajority_reactive_cuff_fraction_r35um_ge2", -1, "Supermajority-reactive cuff within 35 um"),
    ("run4", "round_0003__ca1_pyramidal_reactive_astrocyte_neighbor_fraction", 1, "Literal local reactive-neighbor fraction from April run 4"),
    ("run4", "round_0005__ca1_reactive_local_purity", -1, "Local reactive purity, 60 um"),
    ("run4", "round_0006__ca1_reactive_neuron_exclusion", -1, "Reactive-centered pyramidal-neuron exclusion"),
    ("run5", "round_0002__ca1_pyramidal_near_reactive_astro_fraction_30um", -1, ">=1 reactive astrocyte within 30 um, independent rediscovery"),
    ("run5", "round_0003__ca1_pyramidal_reactive_astro_crowding_fraction_k3_30um", -1, ">=3 reactive astrocytes within 30 um"),
    ("run5", "round_0005__ca1_peripyramidal_astro_reactivity_index_30um", -1, "Peripyramidal astrocyte-lineage reactivity, 30 um"),
    ("run5", "round_0006__ca1_pyramidal_reactive_astro_isolated_fraction_k5_30um", -1, "Reactive-exposed pyramidal neurons with sparse pyramidal neighbors, 30 um"),
]


def load_helpers(path: Path):
    spec = importlib.util.spec_from_file_location("reactive_stats_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_split(input_dir: Path, split: str) -> tuple[pd.DataFrame, dict[tuple[str, str], str]]:
    merged = None
    column_map: dict[tuple[str, str], str] = {}
    metadata = ["donor_id", "slide_name", "slide_path", "slope_zmem0", "max_age_vis", "braak_numeric", "cerad_ordinal", "sex", "sex_binary"]
    for alias, run in RUNS.items():
        frame = pd.read_csv(input_dir / f"{run}_{split}.csv")
        if merged is None:
            merged = frame[metadata].copy()
        for column in frame.columns:
            if column in metadata:
                continue
            renamed = f"{alias}__{column}"
            merged = merged.merge(
                frame[["donor_id", column]].rename(columns={column: renamed}),
                on="donor_id",
                validate="one_to_one",
            )
            column_map[(alias, column)] = renamed
    assert merged is not None
    return merged, column_map


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--helper-script", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-permutations", type=int, default=100_000)
    parser.add_argument("--n-bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument(
        "--measurement-stack",
        default="original April PLIP-based segmentation/classifier surface",
    )
    parser.add_argument("--provenance-note", default="")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    helpers = load_helpers(args.helper_script)

    train, train_map = load_split(args.input_dir, "train")
    test, test_map = load_split(args.input_dir, "test")
    if train_map.keys() != test_map.keys():
        raise RuntimeError("Train/Test feature map mismatch")

    bulk_name = test_map[(BULK_SPEC[0], BULK_SPEC[1])]
    spatial_names = [test_map[(alias, column)] for alias, column, _, _ in SPATIAL_SPECS]
    spatial_signs = np.asarray([row[2] for row in SPATIAL_SPECS], dtype=int)
    labels = {
        test_map[(alias, column)]: label
        for alias, column, _, label in SPATIAL_SPECS
    }

    rows = []
    specs = [(bulk_name, -1, "CA1 reactive-astrocyte fraction", "bulk")]
    specs += [(test_map[(a, c)], s, label, "spatial") for a, c, s, label in SPATIAL_SPECS]
    for split, frame in [("Train", train), ("Test", test)]:
        for offset, (formula, direction, label, role) in enumerate(specs):
            record = helpers.summarize_formula(
                frame,
                formula,
                direction,
                split=split,
                n_bootstrap=args.n_bootstrap,
                seed=args.seed + (0 if split == "Train" else 10_000) + offset,
            )
            record.update({"label": label, "role": role})
            rows.append(record)
    formula_statistics = pd.DataFrame(rows)

    spatial_table, spatial_omnibus = helpers.freedman_lane_matrix(
        test,
        spatial_names,
        spatial_signs,
        n_permutations=args.n_permutations,
        seed=args.seed,
    )
    spatial_table["label"] = spatial_table["formula"].map(labels)
    spatial_table = spatial_table.merge(
        formula_statistics.loc[formula_statistics.split.eq("Test")].drop(
            columns=["split", "label", "role", "partial_r", "declared_direction"]
        ),
        on="formula",
        validate="one_to_one",
    )
    spatial_table.to_csv(args.output_dir / "april_spatial_formula_statistics.csv", index=False)
    formula_statistics.to_csv(args.output_dir / "april_train_test_formula_statistics.csv", index=False)

    bulk_table, bulk_omnibus = helpers.freedman_lane_matrix(
        test,
        [bulk_name],
        np.asarray([-1]),
        n_permutations=args.n_permutations,
        seed=args.seed + 1_000_000,
    )
    p_two = np.asarray(
        [bulk_table.iloc[0].permutation_p_two_sided, spatial_omnibus["maxT_omnibus_p_two_sided"]]
    )
    p_directional = np.asarray(
        [bulk_table.iloc[0].permutation_p_directional, spatial_omnibus["maxT_omnibus_p_directional"]]
    )
    components = pd.DataFrame(
        [
            {
                "component": "Bulk CA1 reactive-astrocyte enrichment",
                "n": len(test),
                "effect_statistic": float(bulk_table.iloc[0].partial_r),
                "effect_statistic_definition": "partial Pearson r",
                "permutation_p_two_sided": p_two[0],
                "permutation_p_directional": p_directional[0],
            },
            {
                "component": "Peripyramidal reactive-astrocyte spatial family",
                "n": len(test),
                "effect_statistic": spatial_omnibus["observed_max_abs_r"],
                "effect_statistic_definition": "maximum absolute partial r across all 14 accepted April spatial entries",
                "permutation_p_two_sided": p_two[1],
                "permutation_p_directional": p_directional[1],
            },
        ]
    )
    components["holm_p_two_sided_across_two_components"] = multipletests(
        p_two, method="holm"
    )[1]
    components["holm_p_directional_across_two_components"] = multipletests(
        p_directional, method="holm"
    )[1]
    components.to_csv(args.output_dir / "april_primary_component_statistics.csv", index=False)

    combined_table, combined_omnibus = helpers.freedman_lane_matrix(
        test,
        [bulk_name, *spatial_names],
        np.asarray([-1, *spatial_signs]),
        n_permutations=args.n_permutations,
        seed=args.seed + 2_000_000,
    )
    combined_table.to_csv(args.output_dir / "april_combined_family_statistics.csv", index=False)

    lodo = []
    for formula, direction, label, role in specs:
        subset = test[["donor_id", formula, helpers.OUTCOME, *helpers.COVARIATES]].dropna().reset_index(drop=True)
        record = helpers.lodo_metrics(subset, formula, direction)
        record.update({"label": label, "role": role})
        lodo.append(record)
    pd.DataFrame(lodo).to_csv(args.output_dir / "april_test_lodo_stability.csv", index=False)

    summary = {
        "measurement_stack": args.measurement_stack,
        "provenance_note": args.provenance_note,
        "n_train": len(train),
        "n_test": len(test),
        "n_spatial_entries": len(spatial_names),
        "settings": {
            "covariates": helpers.COVARIATES,
            "n_permutations": args.n_permutations,
            "n_bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "primary_sidedness": "two-sided",
        },
        "bulk_omnibus": bulk_omnibus,
        "spatial_omnibus": spatial_omnibus,
        "combined_family_omnibus": combined_omnibus,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(components.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Reproduce the TCGA-BRCA train / test split used for the Tubule Score experiments.

Self-contained: reads `tcga_brca_pool.csv` next to this script and writes
`tcga_brca_train.csv` and `tcga_brca_test.csv` next to it. Only the Python
standard library is used. Run from any working directory:

    python split_tcga_brca.py

Input
    tcga_brca_pool.csv   542 TCGA-BRCA cases whose pathology report states a
                         Nottingham tubule score (case_id, tubule_score, note).
                         Available per score: 20 / 118 / 404.

Recipe (seed 42, Python `random`)
    random.seed(42)
    for score in 1, 2, 3:
        cases = sorted(case ids with that score)   # fixed order -> deterministic
        random.shuffle(cases)
        train += cases[:10]
        test  += cases[10:20]

Result: 30 training cases (10 per score) and 30 held-out internal test cases
(10 per score), one slide per patient. The training split is the one every
agent searched on ("Training LOO TCGA-BRCA"); the test split is the one used for
"Internal Validation TCGA-BRCA" and "TCGA-BRCA/Q1.5round*.csv".

With no argument the script also checks the regenerated splits against the
per-case prediction files in the sibling folders and exits non-zero on any
mismatch. Pass --no-check to skip that.
"""
import csv
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 42
N_PER_SCORE = 10


def main():
    groups = defaultdict(list)
    with open(os.path.join(HERE, "tcga_brca_pool.csv"), newline="") as f:
        for row in csv.DictReader(f):
            groups[int(row["tubule_score"])].append(row["case_id"])

    random.seed(SEED)
    train, test = [], []
    for score in sorted(groups):
        cases = sorted(groups[score])
        random.shuffle(cases)
        train += [(c, score) for c in cases[:N_PER_SCORE]]
        test += [(c, score) for c in cases[N_PER_SCORE:2 * N_PER_SCORE]]

    for name, rows in (("tcga_brca_train.csv", train), ("tcga_brca_test.csv", test)):
        rows.sort()
        with open(os.path.join(HERE, name), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["case_id", "tubule_score"])
            w.writerows(rows)
        per = {s: sum(1 for _, t in rows if t == s) for s in sorted(groups)}
        print(f"{name}: n={len(rows)}, per score={per}")

    if "--no-check" in sys.argv:
        return
    parent = os.path.join(HERE, "..")
    refs = [
        ("train", train, f"{parent}/Training LOO TCGA-BRCA/agent_{a}/round20.csv", "case_id")
        for a in range(1, 6)
    ] + [
        ("test", test, f"{parent}/Internal Validation TCGA-BRCA/agent_{a}/round20.csv", "case_id")
        for a in range(1, 6)
    ] + [
        ("test", test, f"{parent}/TCGA-BRCA/Q1.5round{r}.csv", "slide_name") for r in range(1, 6)
    ]
    bad = 0
    for label, rows, path, col in refs:
        if not os.path.exists(path):
            print(f"  [skip] {os.path.relpath(path, HERE)} not found")
            continue
        with open(path, newline="") as f:
            ref = {r[col] for r in csv.DictReader(f)}
        ok = ref == {c for c, _ in rows}
        bad += not ok
        print(f"  [{'OK' if ok else 'MISMATCH'}] {label} == {os.path.relpath(path, HERE)}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

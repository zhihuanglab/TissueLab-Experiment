#!/usr/bin/env python3
"""
Reproduce the HISTAI external cohort (75 cases) used for the Tubule Score
external validation.

Self-contained: reads `histai_breast_pool.csv` next to this script and writes
`histai_batch1.csv`, `histai_batch2.csv`, `histai_batch3.csv` and
`histai_external75.csv` next to it. Only the Python standard library is used.
Run from any working directory:

    python split_histai.py

Input
    histai_breast_pool.csv   858 HISTAI-breast cases whose metadata carries a
                             tubular-differentiation score (case_id, tubule_score).
                             Available per score: 25 / 139 / 694.

Recipe (seed 42, Python `random`)
    The 75 cases were accrued in three batches. Each batch is one stratified
    random draw with the same recipe; the seed is reset once at the start of each
    batch (not once per score), and cases already drawn are removed from the pool:

        random.seed(42)
        for score in 1, 2, 3:
            batch1[score] = random.sample(sorted(pool[score]), 10)
        random.seed(42)
        for score in 1, 2, 3:
            batch2[score] = random.sample(sorted(pool[score] - batch1[score]), 10)
        random.seed(42)
        for score in 1, 2, 3:
            batch3[score] = random.sample(sorted(pool[score] - batch1[score] - batch2[score]), 5)

    `sorted()` fixes the candidate order so the draw is deterministic. Batch 3 is
    5 per score because HISTAI-breast contains only 25 score-1 cases in total and
    20 were already used; 25 / 25 / 25 is the largest balanced cohort the
    collection permits. The union of the three batches is the 75-case external
    cohort ("External Validation histAI"). No HISTAI case was used for training.

With no argument the script also checks the regenerated cohort against the
per-case prediction files in the sibling folder and exits non-zero on any
mismatch. Pass --no-check to skip that.
"""
import csv
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 42
BATCHES = (("histai_batch1.csv", 10), ("histai_batch2.csv", 10), ("histai_batch3.csv", 5))


def write(name, rows):
    rows = sorted(rows)
    with open(os.path.join(HERE, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case_id", "tubule_score"])
        w.writerows(rows)
    per = defaultdict(int)
    for _, s in rows:
        per[s] += 1
    print(f"{name}: n={len(rows)}, per score={dict(sorted(per.items()))}")
    return rows


def main():
    groups = defaultdict(set)
    with open(os.path.join(HERE, "histai_breast_pool.csv"), newline="") as f:
        for row in csv.DictReader(f):
            groups[int(row["tubule_score"])].add(row["case_id"])

    used = defaultdict(set)
    for name, n in BATCHES:
        random.seed(SEED)
        rows = []
        for score in sorted(groups):
            pool = sorted(groups[score] - used[score])
            picked = random.sample(pool, min(n, len(pool)))
            used[score].update(picked)
            rows += [(c, score) for c in picked]
        write(name, rows)
    external = write("histai_external75.csv",
                     [(c, s) for s in groups for c in used[s]])

    if "--no-check" in sys.argv:
        return
    parent = os.path.join(HERE, "..")
    ext_ids = {c for c, _ in external}
    bad = 0
    for a in range(1, 6):
        path = f"{parent}/External Validation histAI/agent_{a}/round20.csv"
        if not os.path.exists(path):
            print(f"  [skip] {os.path.relpath(path, HERE)} not found")
            continue
        with open(path, newline="") as f:
            ref = {r["case_id"] for r in csv.DictReader(f)}
        ok = ref == ext_ids
        bad += not ok
        print(f"  [{'OK' if ok else 'MISMATCH'}] external75 == {os.path.relpath(path, HERE)}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

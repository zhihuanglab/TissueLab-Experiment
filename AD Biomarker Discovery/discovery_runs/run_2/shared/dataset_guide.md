# Dataset Guide: SEA-AD LFB Hippocampus Training Set

## Overview of the dataset

This dataset contains **35 paired whole-slide images (WSIs)** and per-slide **Zarr stores** for LFB-stained hippocampal tissue from the SEA-AD training cohort.

Primary analysis target for biomarker discovery:

- **Outcome:** `slope_zmem0` in `training_cohort.csv`

Key linked data modalities:

1. **Cohort table** (`training_cohort.csv`)
   - donor-level clinical / neuropathology metadata
   - one row per donor / slide
2. **WSI files** (`*.svs`)
   - Aperio pyramidal whole-slide images
3. **Per-slide Zarr stores** (`*.svs.zarr`)
   - cell segmentation geometry
   - 768-d MUSK embeddings
   - cell type classification outputs
   - hippocampal / cortical region annotation polygons

Scale of the cell-level data:

- **35 slides**
- **12,523,007 segmented cells total**
- Per-slide cell counts range from about **163,779** to **528,698** cells
- Core annotated regions present on every slide:
  - `CA1`, `CA2`, `CA3`, `CA4`, `DG`, `EC`, `SB`, `TEC`

---

## File inventory

Top-level `/data` contents:

- `training_cohort.csv` — main donor-level cohort table
- `program.md` — research-program description
- `*.svs` — 35 whole-slide images
- `*.svs.zarr/` — 35 paired Zarr stores, one per slide
- `autoresearch_runs/` — prior run artifacts / logs; not primary raw data
- many `._*` files — macOS AppleDouble sidecar files; **ignore these**

Naming convention:

- donor ID example: `H20.33.001`
- slide / zarr basename example: `H20.33.001-A12-LFB.svs`
- paired zarr example: `H20.33.001-A12-LFB.svs.zarr`

In the cohort table, `slide_name` matches the local zarr directory name exactly.

---

## Cohort table schema (`training_cohort.csv`)

Shape:

- **35 rows × 21 columns**

Columns:

| Column | Type | Notes |
|---|---:|---|
| `donor_id` | string | donor identifier, e.g. `H20.33.001` |
| `slide_name` | string | local zarr basename, e.g. `H20.33.001-A12-LFB.svs.zarr` |
| `slide_path` | string | external/original path, **not** the mounted local `/data` path |
| `slope_zmem0` | float | **primary outcome** |
| `slope_zcasi_irt0` | float | cognitive decline summary |
| `slope_zexf0` | float | executive function decline |
| `slope_zlan0` | float | language decline |
| `slope_zvsp0` | float | visuospatial decline |
| `max_age_vis` | float | confound |
| `MEM_E_last` | float | last-visit memory score |
| `EXF_E_last` | float | last-visit executive score |
| `LAN_E_last` | float | last-visit language score |
| `VSP_E_last` | float | last-visit visuospatial score |
| `braak_label` | string | categorical Braak stage label |
| `cerad_label` | string | categorical CERAD label |
| `overall_ad_neuropath_change` | string | ordered pathology summary |
| `age_at_death` | int | age at death |
| `sex` | string | `Female` / `Male` |
| `cognitive_status` | string | `No dementia` / `Dementia` |
| `braak_numeric` | int | numeric Braak stage |
| `cerad_ordinal` | int | ordinal CERAD |

Observed cohort distributions:

- `sex`: 21 Female, 14 Male
- `cognitive_status`: 19 No dementia, 16 Dementia
- `overall_ad_neuropath_change`: `Not AD`, `Low`, `Intermediate`, `High`
- `braak_numeric` and `cerad_ordinal` already provide analysis-friendly confounds

`sex_binary` is **not present as a column**; if needed for the program's confound set, derive it from `sex`.

### Important cohort notes

- `slide_path` points to the original acquisition location (`/Volumes/...`), not `/data/...`
- For local computation, use:
  - `Path("/data") / slide_name`
- `slide_name` strips to donor ID by removing `-A12-LFB.svs.zarr`

---

## WSI files (`*.svs`)

There are **35 Aperio SVS** slides.

Example properties from one slide:

- dimensions: `(67728, 39040)`
- pyramid levels: `4`
- objective magnification: `20x`
- `aperio.MPP`: `0.5016` or `0.5049` µm/pixel across this cohort

Observed across slides:

- all 35 slides have unique pixel dimensions
- level counts are either **3 or 4**
- MPP values observed: **0.5016** and **0.5049**

Why this matters:

- useful for converting morphology from pixels to approximate microns
- useful if an analysis needs actual image patches from the WSI rather than Zarr-only cell features

---

## Zarr schema (`*.svs.zarr`)

Each slide-level Zarr store has three main groups:

- `SegmentationNode/`
- `ClassificationNode/`
- `CustomAnnotations/`

### 1) `SegmentationNode/`

Core arrays found on every slide:

| Path | Shape | Dtype | Meaning |
|---|---:|---:|---|
| `SegmentationNode/centroids` | `(n_cells, 2)` | `int32` | cell centroids |
| `SegmentationNode/contours` | `(n_cells, 32, 2)` | `int32` | fixed-length contour vertices |
| `SegmentationNode/embedding` | `(n_cells, 768)` | `float16` | MUSK embedding per cell |

Optional / inconsistent:

| Path | Shape | Dtype | Notes |
|---|---:|---:|---|
| `SegmentationNode/probability` | `(n_cells,)` | `float32` | present on **1 slide only** (`H20.33.001-A12-LFB.svs.zarr`) |

Important consistency check:

- for every slide, the leading dimension `n_cells` matches across:
  - centroids
  - contours
  - embeddings
  - class IDs
  - class probabilities

#### Segmentation coordinates

Observed for one representative slide:

- centroid x range: `5 .. 52261`
- centroid y range: `132 .. 37210`

Contours use the same coordinate system as centroids, but with **large absolute coordinates**. Local contour extents are small relative to absolute position (for a sample of 1000 cells, contour vertices were roughly within `±14 px` of the centroid). This is exactly the scenario where area calculations can suffer float32 cancellation.

**Rule:** when computing contour-derived geometry, center each contour first:

```python
c = contour.astype(np.float64)
c = c - c.mean(axis=0)
```

Additional contour observations:

- each contour always has **32 vertices**
- many vertices repeat; effective unique vertex count is often much smaller than 32
- do not assume contours are densely sampled boundaries

---

### 2) `ClassificationNode/`

Common arrays:

| Path | Shape | Dtype | Meaning |
|---|---:|---:|---|
| `ClassificationNode/nuclei_class_id` | `(n_cells,)` | `int32` | predicted class index per cell |
| `ClassificationNode/nuclei_class_probabilities` | `(n_cells, n_classes)` | `float32` | class probability vector |
| `ClassificationNode/nuclei_class_name` | `(n_classes,)` | bytes/string | class names |
| `ClassificationNode/nuclei_class_HEX_color` | `(n_classes,)` | bytes/string | display colors |
| `ClassificationNode/metadata` | scalar | bytes/string | JSON metadata |
| `ClassificationNode/output` | scalar | bytes/string | JSON status summary |

#### Critical class-schema caveat

There are **two class-index schemas** in this cohort:

**Schema A: 34 slides, 9 classes**
```text
0 Negative control
1 Tumor
2 Lymphocyte
3 Pyramidal Neuron
4 Granule Neuron
5 Reactive Astrocyte
6 Astrocyte
7 Oligodendrocyte
8 Corpora Amylacea
```

**Schema B: 1 slide (`H21.33.046-A12-LFB.svs.zarr`), 8 classes**
```text
0 Negative control
1 Pyramidal Neuron
2 Granule Neuron
3 Reactive Astrocyte
4 Astrocyte
5 Oligodendrocyte
6 Corpora Amylacea
7 Lymphocyte
```

Implication:

- **Do not hard-code class IDs globally**
- always map through `nuclei_class_name`
- the `"Tumor"` class exists in 34 slides but had **0 cells globally** in this cohort

Probability arrays are well-formed:

- row sums are ~1.0
- `nuclei_class_id == argmax(nuclei_class_probabilities)` in sampled checks

#### Aggregate class burden across all cells

Global counts across all 35 slides:

- `Oligodendrocyte`: 5,425,491
- `Negative control`: 2,988,444
- `Pyramidal Neuron`: 1,568,643
- `Astrocyte`: 1,367,556
- `Reactive Astrocyte`: 805,363
- `Granule Neuron`: 327,730
- `Lymphocyte`: 26,344
- `Corpora Amylacea`: 13,436
- `Tumor`: 0

Interpretation:

- `Negative control` is a **large fraction** of cells and should usually be excluded or handled explicitly
- `Tumor` is effectively irrelevant here
- biologically meaningful composition features should likely focus on:
  - Pyramidal Neuron
  - Granule Neuron
  - Astrocyte
  - Reactive Astrocyte
  - Oligodendrocyte
  - Lymphocyte
  - Corpora Amylacea

---

### 3) `CustomAnnotations/`

This group contains region polygons. Each annotation is a subgroup like:

- `CA1_01`
- `CA2_01`
- `CA3_01`
- `CA4_01`
- `DG_01`
- `EC_01`
- `SB_01`
- `TEC_01`

Some slides contain extra polygons:

- `EC_02` on 4 slides
- `TEC_02` on 3 slides
- `SB_02` on 1 slide

So the **8 core region bases** are present on all slides, and some regions are split into multiple polygons.

Each annotation subgroup contains:

| Path | Shape | Dtype | Meaning |
|---|---:|---:|---|
| `CustomAnnotations/<region_id>/annotation_json` | scalar | string | W3C-like annotation JSON |
| `CustomAnnotations/<region_id>/comment` | scalar | string | short label, e.g. `CA1-L`, `DG`, `EC` |

Subgroup attributes include:

- `annotation_id`
- `created_at`
- `updated_at`
- `schema`
- `source`
- `wsi_path`

#### Annotation JSON structure

Observed pattern:

- top-level JSON `type` = `"Annotation"`
- polygon lives at:
  - `obj["target"]["selector"]["geometry"]["points"]`
- selector type:
  - `obj["target"]["selector"]["type"] == "POLYGON"`

Point counts vary widely:

- minimum points: **56**
- median: **257**
- maximum: **1006**

#### Critical coordinate-system rule

The annotation polygon coordinates are **16× larger** than the cell coordinate system.

Example from one slide:

- centroid x range: about `0 .. 52k`
- raw annotation x range: about `624k .. 758k`
- after division by 16: about `39k .. 47k`, which matches the cell-coordinate scale

**Rule:** divide annotation coordinates by `16.0` before assigning cells to regions.

```python
pts = np.asarray(obj["target"]["selector"]["geometry"]["points"], dtype=float) / 16.0
```

#### Region-handling recommendation

Treat region membership by **region base name**, not by exact annotation ID:

- combine `EC_01` and `EC_02` into region `EC`
- combine `TEC_01` and `TEC_02` into region `TEC`
- combine `SB_01` and `SB_02` into region `SB`

This will make features portable across slides with split polygons.

---

## Important rules and conventions discovered

### 1) Ignore hidden `._*` files
The dataset contains many macOS AppleDouble sidecar files inside Zarr directories. These can trigger Zarr warnings such as “not recognized as a component of a Zarr hierarchy.”

Practical advice:

- ignore names starting with `._`
- suppress non-fatal warnings if needed during traversal

### 2) Use class names, not numeric class IDs
Because one slide uses a different class ordering and class count, hard-coded IDs are unsafe.

### 3) Divide annotation polygon coordinates by 16
This is necessary for any cell-to-region assignment.

### 4) Center contours before geometry calculations
For area / second-moment / shape features, subtract the contour mean first.

### 5) Expect optional arrays and missing extras
- `SegmentationNode/probability` is not consistently available
- some regions have extra polygons (`_02`), but the 8 core region bases are always present

### 6) Prefer local file paths over `slide_path`
Use `/data/<slide_name>` rather than the external `slide_path` column.

---

## Data structures most relevant to the research objective

For biomarker discovery against `slope_zmem0`, the most useful structures are:

### A. `training_cohort.csv`
Why it matters:

- provides the donor-level outcome
- provides mandatory confounds:
  - `max_age_vis`
  - `braak_numeric`
  - `cerad_ordinal`
  - derived `sex_binary`
- includes useful secondary phenotypes:
  - `overall_ad_neuropath_change`
  - `age_at_death`
  - `cognitive_status`
  - last-visit cognitive domain scores

### B. `SegmentationNode/centroids`
Why it matters:

- enables density, burden, spatial organization, nearest-neighbor, and region occupancy analyses
- base coordinate reference for cell-to-region assignment

### C. `SegmentationNode/contours`
Why it matters:

- direct morphology branch:
  - area
  - perimeter
  - circularity
  - eccentricity / elongation proxies
  - contour roughness

### D. `SegmentationNode/embedding`
Why it matters:

- 768-d MUSK embedding supports:
  - donor-level low-rank summaries
  - per-region embedding averages
  - cluster burden
  - simple linear probes or PCA axes

### E. `ClassificationNode/*`
Why it matters:

- supports composition biomarkers:
  - cell-type fractions
  - burden ratios
  - spatial interaction by type
- use `nuclei_class_name` to map IDs robustly

### F. `CustomAnnotations/*`
Why it matters:

- enables region-specific biomarkers
- region restriction is likely critical biologically, especially for:
  - `CA1`
  - `DG`
  - `EC`
  - hippocampal subfield contrasts

---

## Example code snippets

### Load the cohort table

```python
from pathlib import Path
import pandas as pd

data_root = Path("/data")
cohort = pd.read_csv(data_root / "training_cohort.csv")

# derive program-required sex_binary if needed
cohort["sex_binary"] = (cohort["sex"] == "Male").astype(int)

print(cohort[["donor_id", "slide_name", "slope_zmem0"]].head())
```

### Open a slide Zarr and inspect key arrays

```python
from pathlib import Path
import zarr

slide_name = "H20.33.001-A12-LFB.svs.zarr"
root = zarr.open(Path("/data") / slide_name, mode="r")

centroids = root["SegmentationNode/centroids"]
contours = root["SegmentationNode/contours"]
embedding = root["SegmentationNode/embedding"]
class_ids = root["ClassificationNode/nuclei_class_id"]
class_names = [
    x.decode() if isinstance(x, bytes) else str(x)
    for x in root["ClassificationNode/nuclei_class_name"][:]
]

print(centroids.shape, contours.shape, embedding.shape)
print(class_names)
```

### Map class IDs to class names safely

```python
import numpy as np

class_names = np.array([
    x.decode() if isinstance(x, bytes) else str(x)
    for x in root["ClassificationNode/nuclei_class_name"][:]
])

class_ids = root["ClassificationNode/nuclei_class_id"][:]
cell_type = class_names[class_ids]

# Example: mask biologically relevant cells
mask = np.isin(cell_type, [
    "Pyramidal Neuron",
    "Granule Neuron",
    "Astrocyte",
    "Reactive Astrocyte",
    "Oligodendrocyte",
    "Lymphocyte",
    "Corpora Amylacea",
])
```

### Read region polygons and convert to cell-coordinate scale

```python
import json
import numpy as np

def load_region_polygons(root):
    out = {}
    ann_group = root["CustomAnnotations"]

    for ann_id in ann_group.group_keys():
        raw = ann_group[f"{ann_id}/annotation_json"][()]
        raw = raw.decode() if isinstance(raw, bytes) else raw
        obj = json.loads(raw)

        pts = np.asarray(
            obj["target"]["selector"]["geometry"]["points"],
            dtype=float,
        ) / 16.0  # critical scale correction

        region_base = ann_id.split("_")[0]
        out.setdefault(region_base, []).append(pts)

    return out
```

### Assign cells to a polygon region

```python
import numpy as np
from matplotlib.path import Path as MplPath

xy = root["SegmentationNode/centroids"][:]  # (n_cells, 2)

region_polys = load_region_polygons(root)
ca1_paths = [MplPath(poly) for poly in region_polys["CA1"]]

in_ca1 = np.zeros(len(xy), dtype=bool)
for p in ca1_paths:
    in_ca1 |= p.contains_points(xy)
```

### Compute contour area robustly

```python
import numpy as np

def polygon_area(contour):
    c = contour.astype(np.float64)
    c = c - c.mean(axis=0)  # critical for numerical stability
    x = c[:, 0]
    y = c[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

areas = np.array([polygon_area(c) for c in root["SegmentationNode/contours"][:1000]])
```

### Join donor-level features back to the cohort

```python
from pathlib import Path
import pandas as pd

cohort = pd.read_csv(Path("/data") / "training_cohort.csv")
cohort["sex_binary"] = (cohort["sex"] == "Male").astype(int)

# suppose donor_scores is a DataFrame with columns: donor_id, biomarker
df = cohort.merge(donor_scores, on="donor_id", how="left")

analysis_cols = [
    "biomarker",
    "slope_zmem0",
    "max_age_vis",
    "braak_numeric",
    "cerad_ordinal",
    "sex_binary",
]
analysis_df = df[analysis_cols]
```

---

## Recommended analysis defaults

For portable biomarker engineering, a good default preprocessing policy is:

1. use `training_cohort.csv` as the donor roster
2. open each slide via `/data/<slide_name>`
3. load class names dynamically for each slide
4. exclude or explicitly model `Negative control`
5. ignore `Tumor`
6. combine multiple polygons per region base
7. divide annotation coordinates by 16 before region assignment
8. center contours before shape calculations
9. return `NaN` / `None` gracefully if a region is missing or has too few cells

---

## Bottom line

The strongest, safest raw ingredients for this project are:

- **cohort outcome/confounds** from `training_cohort.csv`
- **region-restricted cell composition** from `ClassificationNode`
- **morphology features** from `SegmentationNode/contours`
- **spatial organization** from `SegmentationNode/centroids`
- **simple embedding summaries** from `SegmentationNode/embedding`

The two biggest implementation pitfalls are:

1. **annotation polygons must be divided by 16**
2. **cell class IDs are not globally stable across slides**

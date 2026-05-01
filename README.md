# Phlow

Phlow is a small Python toolkit for organizing and analyzing flow cytometry `.fcs`
experiments. It can sort raw experiment files into label folders, read gated GFP
and mCherry/RFP data with FlowCal, generate per-strain summary plots, and create
comparison plots across strains.

The project is built around simple scripts, so it works both as a command-line
workflow and as importable Python utilities.

## Features

- Organize raw `.fcs` files into one subfolder per strain/label.
- Normalize `.fcs` numbering inside label folders.
- Optionally reverse file numbering when explicitly requested.
- Validate expected file counts before processing.
- Support standard condition layouts and triplicate experiments.
- Read `.fcs` files with FlowCal.
- Gate events using forward/side scatter and GFP/RFP channel thresholds.
- Compile GFP, mCherry/RFP, and GFP/RFP-normalized data.
- Save seaborn histogram/KDE plots for each label.
- Save GFP and mCherry scatter plots versus light input.
- Compare multiple labels on shared GFP and mCherry scatter plots.
- Regenerate customized plots for a single already-processed label.
- Use custom display labels without changing folder/strain names.
- Configure scatter plot y-axis limits from the CLI.
- Estimate mutual information/channel capacity from binned flow data.

## Repository Layout

```text
.
├── flow_object.py      # Core flow_unit class, plotting, metrics, capacity code
├── flow_organize.py    # File organization and renaming helpers
├── flow_pipeline.py    # Main organization + per-label output pipeline
├── flow_compare.py     # Multi-label comparison and individualized replotting
└── README.md
```

## Dependencies

Phlow expects a Python environment with:

```text
numpy
matplotlib
seaborn
FlowCal
```

Install dependencies with your preferred environment manager. For example:

```bash
pip install numpy matplotlib seaborn FlowCal
```

## Data Assumptions

Phlow expects `.fcs` filenames to end with a numeric suffix:

```text
sample-name-1.fcs
sample-name-2.fcs
sample-name-3.fcs
...
```

For `num_cond=4` without triplicates, each label has four `.fcs` files.

For `num_cond=4` with triplicates, each label has twelve `.fcs` files: three
replicates for each condition.

By default, the four light input labels are:

```text
0, 21, 52, 208 µW/cm^2
```

Scatter plots use evenly spaced categorical x positions so the `208` condition
does not visually sit far away from the other conditions. The tick labels still
show the real light input values.

## Main Pipeline

Use `flow_pipeline.py` when you want to organize raw files and generate all
basic outputs.

### Raw `.fcs` Files in Root Folder

If your root folder contains raw `.fcs` files that still need to be organized:

```bash
python flow_pipeline.py /path/to/experiment \
  --labels strainA strainB strainC
```

This will:

1. Validate the expected `.fcs` file count.
2. Create one folder per label.
3. Move files into label folders.
4. Normalize numbering inside each label folder.
5. Read and compile flow data.
6. Save plots inside each label folder.

### Already Organized Label Folders

If your root folder already contains label folders, labels are optional:

```bash
python flow_pipeline.py /path/to/experiment
```

Only subfolders containing `.fcs` files are treated as labels.

### Triplicates

Use `--triplicate` when each condition has three replicate `.fcs` files:

```bash
python flow_pipeline.py /path/to/experiment \
  --labels strainA strainB \
  --triplicate
```

With `num_cond=4`, this expects `12` files per label.

### Custom Condition Count and Light Inputs

```bash
python flow_pipeline.py /path/to/experiment \
  --labels strainA strainB \
  --num-cond 5 \
  --light-inputs 0,10,25,50,100
```

The number of light inputs must match `--num-cond`.

### Custom Plot Labels

Folder names are used for file loading. Plot labels are used only in legends and
display text.

```bash
python flow_pipeline.py /path/to/experiment \
  --labels strain_folder_A strain_folder_B \
  --plot-labels "WT" "Mutant"
```

If `--plot-labels` is omitted, Phlow uses the folder/strain labels.

### Scatter Y-Limits

Defaults:

- GFP: `3,4.5`
- mCherry: `0,2`

Override them with:

```bash
python flow_pipeline.py /path/to/experiment \
  --labels strainA strainB \
  --gfp-ylim 2.8,4.8 \
  --mcherry-ylim 0,2.5
```

### Optional Reverse Numbering

File numbering is not reversed unless you explicitly request it:

```bash
python flow_pipeline.py /path/to/experiment \
  --labels strainA strainB \
  --reverse-numbering
```

## Per-Label Outputs

For each label folder, the pipeline saves:

```text
histograms.svg
gfp_vs_light_input.svg
mcherry_vs_light_input.svg
```

The histogram file contains three KDE panels:

- mCherry/RFP
- GFP
- normalized GFP/RFP

Scatter plots are saved as SVG at `dpi=300`.

## Compare Labels

Use `flow_compare.py` when your experiment is already organized into label
folders and you want comparison plots.

```bash
python flow_compare.py /path/to/experiment \
  --labels strainA strainB strainC
```

This saves root-level comparison plots:

```text
compare_gfp_vs_light_input_strainA__strainB__strainC.svg
compare_mcherry_vs_light_input_strainA__strainB__strainC.svg
```

Comparison filenames include the compared labels to avoid overwriting previous
comparisons.

In comparison plots:

- GFP points are green.
- mCherry points are red.
- Labels are distinguished by marker shape.
- Legends use `--plot-labels` when provided.

Example with display labels:

```bash
python flow_compare.py /path/to/experiment \
  --labels strain_folder_A strain_folder_B \
  --plot-labels "WT" "Mutant"
```

## Individualized Replotting

`flow_compare.py` also works as an individualized plot processor. If you pass
exactly one label, it regenerates that label's histograms and scatter plots in
the label folder instead of writing comparison plots to the root folder.

```bash
python flow_compare.py /path/to/experiment \
  --labels strainA \
  --gfp-ylim 3.1,4.2 \
  --mcherry-ylim 0.2,1.7
```

This is useful when the global pipeline settings were not ideal for one strain
and you want to adjust only that strain's plots.

## Python API

You can call the pipeline directly from Python:

```python
from flow_pipeline import run_flow_experiment

outputs = run_flow_experiment(
    root_folder="/path/to/experiment",
    labels=["strainA", "strainB"],
    num_cond=4,
    triplicate=False,
    gain=8,
    light_inputs=[0, 21, 52, 208],
    gfp_ylim=(3, 4.5),
    mcherry_ylim=(0, 2),
    plot_labels=["WT", "Mutant"],
)
```

Compare labels from Python:

```python
from flow_compare import compare_labels

outputs = compare_labels(
    root_folder="/path/to/experiment",
    labels=["strainA", "strainB"],
    plot_labels=["WT", "Mutant"],
)
```

## Core `flow_unit`

`flow_object.py` defines the lower-level `flow_unit` object for one strain or
experiment label.

It supports:

- setting `.fcs` file paths
- setting gain values
- configuring number of conditions
- triplicate merging
- log or linear data modes
- reading gated GFP and RFP data
- gain-corrected GFP compilation
- normalized GFP/RFP data
- population metrics
- histograms and scatter plots
- channel capacity calculations

Basic usage:

```python
from flow_object import flow_unit

unit = flow_unit("strainA", "WT")
unit.set_paths("/path/to/experiment/strainA/sample-")
unit.set_gain(8)
unit.read_data()
unit.compile_data()

rfp_metrics, gfp_metrics, corrected_metrics = unit.compute_pop_metrics()
capacity = unit.compute_CC()
```

## File Organization Helpers

`flow_organize.py` provides reusable helpers:

```python
from flow_organize import (
    organize_fcs_files,
    normalize_filenames,
    normalize_in_all_subfolders,
    reverse_fcs_numbering,
)
```

Available helpers:

- `organize_fcs_files(...)`: move files into label folders by trailing number.
- `normalize_filenames(...)`: rename files in one folder to sequential suffixes.
- `normalize_in_all_subfolders(...)`: normalize every label folder.
- `reverse_fcs_numbering(...)`: reverse suffix numbering using temporary UUID
  names to avoid collisions.

## Notes

- Plot titles are intentionally omitted from generated SVGs.
- Scatter plot x-positions are categorical, but tick labels show true light
  input values.
- The default gain is `8`, matching the existing GFP correction formula
  `GFP * 8 / gain`.
- `FlowCal` must be installed and able to read your `.fcs` files.


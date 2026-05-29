# Phlow

Phlow is a small Python toolkit for organizing and analyzing flow cytometry `.fcs`
experiments. It can sort raw experiment files into label folders, read gated GFP
and mCherry/RFP data with FlowCal, generate per-strain summary plots, and create
comparison plots across strains.

The project is now packaged as `Phlow`. After installation, the top-level
terminal command is `phlow` or `Phlow`, with subcommands that call the original
pipeline and comparison code.

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
├── pyproject.toml
├── phlow/
│   ├── cli.py           # Top-level phlow / Phlow command-line interface
│   ├── flow_object.py   # Core flow_unit class, plotting, metrics, capacity code
│   ├── flow_organize.py # File organization and renaming helpers
│   ├── flow_pipeline.py # Main organization + per-label output pipeline
│   └── flow_compare.py  # Multi-label comparison and individualized replotting
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

Install the package in editable mode from the repository root:

```bash
pip install -e .
```

This installs two equivalent command names:

```bash
phlow --help
Phlow --help
```

You can also run the package without installing the console script:

```bash
python -m phlow --help
```

## Command-Line Layer

The package CLI lives in `phlow/cli.py` and provides three subcommands:

```text
phlow run
phlow compare
phlow g-compare
```

`phlow run` and `phlow compare` accept `--address`, which is the path to one
experiment's mother folder. If `--address` is omitted, Phlow uses the current
working directory. `phlow g-compare` accepts `--addresses`, a list of two or
more experiment mother folders.

The commands share the main analysis options:

```text
--labels label1 label2 ...
--num-cond 4
--num_cond 4
--triplicate
--gain 8
--light-inputs 0,21,52,208
--plot-labels "WT" "Mutant"
--gfp-ylim 2.8,4.8
--mcherry-ylim 0,2.5
```

`--num-cond` and `--num_cond` are equivalent.

The CLI is a thin layer over the existing modules:

- `phlow run` calls `phlow.flow_pipeline.run_flow_experiment(...)`.
- `phlow compare` calls `phlow.flow_compare.compare_labels(...)`.
- `phlow g-compare` calls `phlow.global_compare.global_compare(...)`.
- `phlow run --reverse-numbering` also uses the file-renaming helpers in
  `phlow.flow_organize`.

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

Use `phlow run` when you want to organize raw files and generate all basic
per-label outputs.

### Raw `.fcs` Files in Root Folder

If your root folder contains raw `.fcs` files that still need to be organized:

```bash
phlow run --address /path/to/experiment \
  --labels strainA strainB strainC
```

This will:

1. Validate the expected `.fcs` file count.
2. Create one folder per label.
3. Move files into label folders.
4. Normalize numbering inside each label folder.
5. Read and compile flow data.
6. Save plots inside each label folder.

If `--address` is not provided, the same command runs against the current
working directory:

```bash
cd /path/to/experiment
phlow run --labels strainA strainB strainC
```

### Already Organized Label Folders

If your root folder already contains label folders, labels are optional:

```bash
phlow run --address /path/to/experiment
```

Only subfolders containing `.fcs` files are treated as labels. In this mode,
Phlow skips organization and normalization and only regenerates plots for each
label folder.

### Triplicates

Use `--triplicate` when each condition has three replicate `.fcs` files:

```bash
phlow run --address /path/to/experiment \
  --labels strainA strainB \
  --triplicate
```

With `num_cond=4`, this expects `12` files per label.

### Custom Condition Count and Light Inputs

```bash
phlow run --address /path/to/experiment \
  --labels strainA strainB \
  --num-cond 5 \
  --light-inputs 0,10,25,50,100
```

The number of light inputs must match `--num-cond`.

### Custom Plot Labels

Folder names are used for file loading. Plot labels are used only in legends and
display text.

```bash
phlow run --address /path/to/experiment \
  --labels strain_folder_A strain_folder_B \
  --plot-labels "WT" "Mutant"
```

If `--plot-labels` is omitted, Phlow uses the folder/strain labels.

### Scatter Y-Limits

By default, Phlow calculates scatter y-limits from the plotted data and adds a
small clearance above and below the observed values.

You can override the automatic limits with:

```bash
phlow run --address /path/to/experiment \
  --labels strainA strainB \
  --gfp-ylim 2.8,4.8 \
  --mcherry-ylim 0,2.5
```

### Optional Reverse Numbering

File numbering is not reversed unless you explicitly request it:

```bash
phlow run --address /path/to/experiment \
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

Use `phlow compare` when your experiment is already organized into label
folders and you want comparison plots.

```bash
phlow compare --address /path/to/experiment \
  --labels strainA strainB strainC
```

With no `--address`, comparison also defaults to the current working directory:

```bash
cd /path/to/experiment
phlow compare --labels strainA strainB strainC
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
phlow compare --address /path/to/experiment \
  --labels strain_folder_A strain_folder_B \
  --plot-labels "WT" "Mutant"
```

## Individualized Replotting

`phlow compare` also works as an individualized plot processor. If you pass
exactly one label, it regenerates that label's histograms and scatter plots in
the label folder instead of writing comparison plots to the root folder.

```bash
phlow compare --address /path/to/experiment \
  --labels strainA \
  --gfp-ylim 3.1,4.2 \
  --mcherry-ylim 0.2,1.7
```

This is useful when the global pipeline settings were not ideal for one strain
and you want to adjust only that strain's plots.

## Global Comparison Across Days

Use `phlow g-compare` when you want to compare flow data collected in different
experiment folders, such as measurements from different days.

The command expects one label per address, in the same order:

```bash
phlow g-compare \
  --addresses /path/to/day1 /path/to/day2 /path/to/day3 \
  --labels strainA strainA strainA
```

This loads:

```text
/path/to/day1/strainA
/path/to/day2/strainA
/path/to/day3/strainA
```

The command saves identical global comparison outputs into a comparison-specific
subfolder inside `Global comparisons` for every address. The subfolder name is
built from the ordered labels, such as `strainA-strainA-strainA`:

```text
/path/to/day1/Global comparisons/strainA-strainA-strainA/
/path/to/day2/Global comparisons/strainA-strainA-strainA/
/path/to/day3/Global comparisons/strainA-strainA-strainA/
```

Global comparison outputs include:

```text
global_compare_gfp_vs_light_input_*.svg
global_compare_mcherry_vs_light_input_*.svg
global_compare_gfp_histograms_*.svg
```

### Global Gains and Triplicates

By default, every address uses gain `8` and `triplicate=False`.

To provide per-address gains:

```bash
phlow g-compare \
  --addresses /path/to/day1 /path/to/day2 /path/to/day3 \
  --labels strainA strainA strainA \
  --gain 1 8 64
```

To provide per-address triplicate settings:

```bash
phlow g-compare \
  --addresses /path/to/day1 /path/to/day2 /path/to/day3 \
  --labels strainA strainA strainA \
  --triplicates true false true
```

The number of values passed to `--labels`, `--gain`, and `--triplicates` must
match the number of values passed to `--addresses`.

### Global Plot Labels

If `--plot-labels` is omitted, global comparison legends use the label plus the
experiment folder name, for example `strainA (day1)`.

You can override those display labels:

```bash
phlow g-compare \
  --addresses /path/to/day1 /path/to/day2 \
  --labels strainA strainA \
  --plot-labels "Day 1" "Day 2"
```

You can also customize the scatter plot y-axis labels:

```bash
phlow g-compare \
  --addresses /path/to/day1 /path/to/day2 \
  --labels strainA strainA \
  --GFP-y-label "Mean GFP" \
  --mCherry-y-label "Mean mCherry"
```

## Python API

You can call the pipeline directly from Python:

```python
from phlow.flow_pipeline import run_flow_experiment

outputs = run_flow_experiment(
    root_folder="/path/to/experiment",
    labels=["strainA", "strainB"],
    num_cond=4,
    triplicate=False,
    gain=8,
    light_inputs=[0, 21, 52, 208],
    plot_labels=["WT", "Mutant"],
)
```

Compare labels from Python:

```python
from phlow.flow_compare import compare_labels

outputs = compare_labels(
    root_folder="/path/to/experiment",
    labels=["strainA", "strainB"],
    plot_labels=["WT", "Mutant"],
)
```

Global comparison from Python:

```python
from phlow.global_compare import global_compare

outputs = global_compare(
    addresses=["/path/to/day1", "/path/to/day2"],
    labels=["strainA", "strainA"],
    gains=[8, 64],
    triplicates=[False, True],
    plot_labels=["Day 1", "Day 2"],
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
from phlow.flow_object import flow_unit

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
from phlow.flow_organize import (
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

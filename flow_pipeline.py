import argparse
import os
import re

from flow_organize import (
    normalize_in_all_subfolders,
    organize_fcs_files,
    reverse_fcs_numbering,
)


DEFAULT_LIGHT_INPUTS = [0, 21, 52, 208]
DEFAULT_GFP_YLIM = (3, 4.5)
DEFAULT_MCHERRY_YLIM = (0, 2)


def _count_fcs_files(folder):
    return sum(1 for name in os.listdir(folder) if name.lower().endswith(".fcs"))


def _label_folder(root_folder, label):
    return os.path.join(root_folder, str(label))


def _fcs_subfolders(root_folder):
    return [
        name
        for name in sorted(os.listdir(root_folder))
        if os.path.isdir(os.path.join(root_folder, name))
        and _count_fcs_files(os.path.join(root_folder, name)) > 0
    ]


def _infer_fcs_prefix(folder):
    pattern = re.compile(r"(.+-)[0-9]+\.fcs$", re.IGNORECASE)
    for filename in sorted(os.listdir(folder)):
        match = pattern.match(filename)
        if match:
            return os.path.join(folder, match.group(1))
    raise ValueError(f"No normalized .fcs files ending in '-#.fcs' found in {folder}")


def _condition_labels(num_cond, light_inputs=None):
    if light_inputs is None:
        light_inputs = DEFAULT_LIGHT_INPUTS if num_cond == 4 else list(range(1, num_cond + 1))
    if len(light_inputs) != num_cond:
        raise ValueError("light_inputs must have one value per condition")
    return [str(value) for value in light_inputs]


def _resolve_plot_labels(labels, plot_labels=None):
    if plot_labels is None:
        return {str(label): str(label) for label in labels}
    if len(plot_labels) != len(labels):
        raise ValueError("plot_labels must have the same number of entries as labels")
    return {
        str(label): str(plot_label)
        for label, plot_label in zip(labels, plot_labels)
    }


def _categorical_x_positions(num_points):
    return list(range(1, (num_points * 2) + 1, 2))


def _save_histograms(unit, label_folder, label):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(6.5, 8.0), constrained_layout=True)
    unit.plot_histogram(axes=axes, is_legend=True, is_norm=True)

    axes[0].set_xlabel("log10(RFP)")
    axes[0].set_ylabel("Density")
    axes[1].set_xlabel("log10(GFP)")
    axes[1].set_ylabel("Density")
    axes[2].set_xlabel("log10(GFP/RFP)")
    axes[2].set_ylabel("Density")

    output_path = os.path.join(label_folder, "histograms.svg")
    fig.savefig(output_path, format="svg", dpi=300)
    plt.close(fig)
    return output_path


def _save_channel_light_scatter(
    light_inputs,
    means,
    output_path,
    ylabel,
    color,
    ylim,
):
    import matplotlib.pyplot as plt

    x_positions = _categorical_x_positions(len(light_inputs))
    fig, ax = plt.subplots(figsize=(5.5, 4.0), constrained_layout=True)
    ax.scatter(x_positions, means, color=color, marker="o")
    ax.set_xlabel(r"Green light intensity ($\mu$W/cm$^2$)")
    ax.set_ylabel(ylabel)
    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(value) for value in light_inputs])
    ax.set_ylim(ylim)

    fig.savefig(output_path, format="svg", dpi=300)
    plt.close(fig)
    return output_path


def _save_light_scatter_outputs(
    unit,
    label_folder,
    label,
    light_inputs,
    gfp_ylim=DEFAULT_GFP_YLIM,
    mcherry_ylim=DEFAULT_MCHERRY_YLIM,
):
    rfp_metrics, gfp_metrics, _ = unit.compute_pop_metrics()
    return {
        "gfp_vs_light_input": _save_channel_light_scatter(
            light_inputs=light_inputs,
            means=gfp_metrics["Mean"],
            output_path=os.path.join(label_folder, "gfp_vs_light_input.svg"),
            ylabel="Mean log10(GFP)",
            color="green",
            ylim=gfp_ylim,
        ),
        "mcherry_vs_light_input": _save_channel_light_scatter(
            light_inputs=light_inputs,
            means=rfp_metrics["Mean"],
            output_path=os.path.join(label_folder, "mcherry_vs_light_input.svg"),
            ylabel="Mean log10(mCherry)",
            color="red",
            ylim=mcherry_ylim,
        ),
    }


def load_label_flow_unit(
    label,
    root_folder,
    num_cond,
    triplicate,
    gain,
    cond_labels,
    plot_label=None,
):
    from flow_object import flow_unit

    folder = _label_folder(root_folder, label)
    unit = flow_unit(str(label), str(plot_label if plot_label is not None else label))
    unit.set_number_of_conditions(num_cond)
    unit.cond_labels = cond_labels
    unit.set_triplicates(triplicate)
    unit.set_paths(_infer_fcs_prefix(folder))
    unit.set_gain(gain)
    unit.read_data()
    unit.compile_data()
    return unit


def run_flow_experiment(
    root_folder,
    labels=None,
    num_cond=4,
    triplicate=False,
    reverse_numbering=False,
    gain=8,
    light_inputs=None,
    gfp_ylim=DEFAULT_GFP_YLIM,
    mcherry_ylim=DEFAULT_MCHERRY_YLIM,
    plot_labels=None,
):
    """
    Organize a flow cytometry experiment and save basic plots for each label.

    Parameters
    ----------
    root_folder : str
        Folder containing the raw .fcs files.
    labels : list[str], optional
        Folder/condition labels. Required when raw .fcs files still need to be
        organized from root_folder. If omitted, existing subfolders in
        root_folder are treated as labels and only outputs are generated.
    num_cond : int, optional
        Number of experimental conditions per label.
    triplicate : bool, optional
        If True, expects three .fcs files per condition.
    reverse_numbering : bool, optional
        If True, reverses .fcs numbering inside each label folder after
        organizing and normalizing.
    gain : int or float, optional
        Flow cytometer gain passed to flow_unit. Defaults to 8, which leaves
        GFP gain scaling unchanged in the existing correction formula.
    light_inputs : list[int | float], optional
        Light input values used on plots. Defaults to [0, 21, 52, 208] when
        num_cond is 4, otherwise [1, ..., num_cond].
    gfp_ylim : tuple[float, float], optional
        Y-axis limits for GFP scatter plots. Default: (3, 4.5).
    mcherry_ylim : tuple[float, float], optional
        Y-axis limits for mCherry scatter plots. Default: (0, 2).
    plot_labels : list[str], optional
        Labels to use in plot titles instead of folder/strain labels.

    Returns
    -------
    dict
        Mapping from label to generated output file paths.
    """
    if not os.path.isdir(root_folder):
        raise ValueError(f"root_folder does not exist or is not a directory: {root_folder}")

    existing_label_folders = _fcs_subfolders(root_folder)
    if labels is None:
        if not existing_label_folders:
            raise ValueError(
                "labels are required when root_folder does not already contain label subfolders"
            )
        labels = existing_label_folders
    else:
        labels = [str(label) for label in labels]
        if not labels:
            raise ValueError("labels must contain at least one label")

    labels_are_already_organized = all(
        os.path.isdir(_label_folder(root_folder, label)) for label in labels
    )

    if light_inputs is None:
        light_inputs = DEFAULT_LIGHT_INPUTS if num_cond == 4 else list(range(1, num_cond + 1))
    else:
        light_inputs = list(light_inputs)
    cond_labels = _condition_labels(num_cond, light_inputs)
    plot_label_map = _resolve_plot_labels(labels, plot_labels)
    files_per_label = num_cond * (3 if triplicate else 1)

    if labels_are_already_organized:
        for label in labels:
            folder = _label_folder(root_folder, label)
            actual_in_label = _count_fcs_files(folder)
            if actual_in_label != files_per_label:
                raise ValueError(
                    f"Expected {files_per_label} .fcs files in {folder}, found {actual_in_label}"
                )
    else:
        expected_total = files_per_label * len(labels)
        actual_total = _count_fcs_files(root_folder)

        if actual_total != expected_total:
            raise ValueError(
                f"Expected {expected_total} .fcs files in {root_folder} "
                f"({files_per_label} per label), found {actual_total}"
            )

        organize_fcs_files(root_folder, labels, num_conds=files_per_label)
        normalize_in_all_subfolders(root_folder)

    if reverse_numbering:
        for label in labels:
            reverse_fcs_numbering(_label_folder(root_folder, label))

    outputs = {}
    for label in labels:
        folder = _label_folder(root_folder, label)
        expected_in_label = files_per_label
        actual_in_label = _count_fcs_files(folder)
        if actual_in_label != expected_in_label:
            raise ValueError(
                f"Expected {expected_in_label} .fcs files in {folder}, found {actual_in_label}"
            )

        plot_label = plot_label_map[str(label)]
        unit = load_label_flow_unit(
            label,
            root_folder,
            num_cond,
            triplicate,
            gain,
            cond_labels,
            plot_label=plot_label,
        )

        outputs[str(label)] = {
            "histograms": _save_histograms(unit, folder, plot_label),
        }
        outputs[str(label)].update(
            _save_light_scatter_outputs(
                unit,
                folder,
                plot_label,
                light_inputs,
                gfp_ylim=gfp_ylim,
                mcherry_ylim=mcherry_ylim,
            )
        )

    return outputs


def _parse_light_inputs(value):
    try:
        return [float(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "light inputs must be a comma-separated list of numbers"
        ) from exc


def _parse_ylim(value):
    try:
        parts = [float(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "ylim must be two comma-separated numbers, e.g. '3,4.5'"
        ) from exc
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            "ylim must be two comma-separated numbers, e.g. '3,4.5'"
        )
    if parts[0] >= parts[1]:
        raise argparse.ArgumentTypeError("ylim minimum must be smaller than maximum")
    return tuple(parts)


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Organize flow cytometry .fcs files and generate basic SVG outputs."
    )
    parser.add_argument(
        "root_folder",
        help="Root folder containing raw .fcs files or existing label subfolders.",
    )
    parser.add_argument(
        "-l",
        "--labels",
        nargs="+",
        help=(
            "Label names. Required when organizing raw .fcs files from root_folder. "
            "Optional when root_folder already contains label subfolders."
        ),
    )
    parser.add_argument(
        "-n",
        "--num-cond",
        type=int,
        default=4,
        help="Number of conditions per label. Default: 4.",
    )
    parser.add_argument(
        "-t",
        "--triplicate",
        action="store_true",
        help="Expect three .fcs files per condition.",
    )
    parser.add_argument(
        "-r",
        "--reverse-numbering",
        action="store_true",
        help="Reverse .fcs numbering inside each label folder after normalization.",
    )
    parser.add_argument(
        "-g",
        "--gain",
        type=float,
        default=8,
        help="Flow cytometer gain passed to flow_unit. Default: 8.",
    )
    parser.add_argument(
        "--light-inputs",
        type=_parse_light_inputs,
        help="Comma-separated light input values, e.g. '0,21,52,208'.",
    )
    parser.add_argument(
        "--plot-labels",
        nargs="+",
        help=(
            "Labels to show in plots. Must match --labels count when provided. "
            "Defaults to the folder/strain labels."
        ),
    )
    parser.add_argument(
        "--gfp-ylim",
        type=_parse_ylim,
        default=DEFAULT_GFP_YLIM,
        help="GFP scatter y-axis limits as min,max. Default: 3,4.5.",
    )
    parser.add_argument(
        "--mcherry-ylim",
        type=_parse_ylim,
        default=DEFAULT_MCHERRY_YLIM,
        help="mCherry scatter y-axis limits as min,max. Default: 0,2.",
    )
    return parser


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    outputs = run_flow_experiment(
        root_folder=args.root_folder,
        labels=args.labels,
        num_cond=args.num_cond,
        triplicate=args.triplicate,
        reverse_numbering=args.reverse_numbering,
        gain=args.gain,
        light_inputs=args.light_inputs,
        gfp_ylim=args.gfp_ylim,
        mcherry_ylim=args.mcherry_ylim,
        plot_labels=args.plot_labels,
    )

    for label, paths in outputs.items():
        print(f"{label}:")
        for output_name, output_path in paths.items():
            print(f"  {output_name}: {output_path}")


if __name__ == "__main__":
    main()


__all__ = [
    "DEFAULT_GFP_YLIM",
    "DEFAULT_MCHERRY_YLIM",
    "_resolve_plot_labels",
    "load_label_flow_unit",
    "run_flow_experiment",
]

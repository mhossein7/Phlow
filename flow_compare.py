import argparse
import os
import re

from flow_pipeline import (
    DEFAULT_LIGHT_INPUTS,
    DEFAULT_GFP_YLIM,
    DEFAULT_MCHERRY_YLIM,
    _categorical_x_positions,
    _condition_labels,
    _count_fcs_files,
    _label_folder,
    _parse_light_inputs,
    _parse_ylim,
    _resolve_plot_labels,
    load_label_flow_unit,
    run_flow_experiment,
)


def compare_labels(
    root_folder,
    labels,
    num_cond=4,
    triplicate=False,
    gain=8,
    light_inputs=None,
    gfp_ylim=DEFAULT_GFP_YLIM,
    mcherry_ylim=DEFAULT_MCHERRY_YLIM,
    plot_labels=None,
):
    """
    Save individualized plots for one label, or root-level comparison plots for
    multiple labels.
    """
    if not os.path.isdir(root_folder):
        raise ValueError(f"root_folder does not exist or is not a directory: {root_folder}")
    if not labels:
        raise ValueError("labels must contain at least one label")

    labels = [str(label) for label in labels]
    if light_inputs is None:
        light_inputs = DEFAULT_LIGHT_INPUTS if num_cond == 4 else list(range(1, num_cond + 1))
    else:
        light_inputs = list(light_inputs)

    if len(labels) == 1:
        return run_flow_experiment(
            root_folder=root_folder,
            labels=labels,
            num_cond=num_cond,
            triplicate=triplicate,
            gain=gain,
            light_inputs=light_inputs,
            gfp_ylim=gfp_ylim,
            mcherry_ylim=mcherry_ylim,
            plot_labels=plot_labels,
        )

    cond_labels = _condition_labels(num_cond, light_inputs)
    plot_label_map = _resolve_plot_labels(labels, plot_labels)
    files_per_label = num_cond * (3 if triplicate else 1)

    units = {}
    for label in labels:
        folder = _label_folder(root_folder, label)
        if not os.path.isdir(folder):
            raise ValueError(f"Label folder does not exist: {folder}")
        actual_files = _count_fcs_files(folder)
        if actual_files != files_per_label:
            raise ValueError(
                f"Expected {files_per_label} .fcs files in {folder}, found {actual_files}"
            )
        units[label] = load_label_flow_unit(
            label=label,
            root_folder=root_folder,
            num_cond=num_cond,
            triplicate=triplicate,
            gain=gain,
            cond_labels=cond_labels,
            plot_label=plot_label_map[label],
        )

    outputs = {
        "gfp_comparison": _save_comparison_scatter(
            root_folder=root_folder,
            units=units,
            plot_label_map=plot_label_map,
            light_inputs=light_inputs,
            metric_name="GFP",
            ylabel="Mean log10(GFP)",
            ylim=gfp_ylim,
            output_filename=f"compare_gfp_vs_light_input_{_labels_filename_suffix(labels)}.svg",
        ),
        "mcherry_comparison": _save_comparison_scatter(
            root_folder=root_folder,
            units=units,
            plot_label_map=plot_label_map,
            light_inputs=light_inputs,
            metric_name="mCherry",
            ylabel="Mean log10(mCherry)",
            ylim=mcherry_ylim,
            output_filename=f"compare_mcherry_vs_light_input_{_labels_filename_suffix(labels)}.svg",
        ),
    }
    return outputs


def _labels_filename_suffix(labels):
    safe_labels = []
    for label in labels:
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(label)).strip("_")
        safe_labels.append(safe or "label")
    return "__".join(safe_labels)


def _save_comparison_scatter(
    root_folder,
    units,
    plot_label_map,
    light_inputs,
    metric_name,
    ylabel,
    ylim,
    output_filename,
):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.0, 4.2), constrained_layout=True)
    x_positions = _categorical_x_positions(len(light_inputs))
    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]
    color = "green" if metric_name == "GFP" else "red"

    for idx, (label, unit) in enumerate(units.items()):
        rfp_metrics, gfp_metrics, _ = unit.compute_pop_metrics()
        means = gfp_metrics["Mean"] if metric_name == "GFP" else rfp_metrics["Mean"]
        ax.scatter(
            x_positions,
            means,
            color=color,
            marker=markers[idx % len(markers)],
            label=plot_label_map[label],
        )

    ax.set_xlabel(r"Green light intensity ($\mu$W/cm$^2$)")
    ax.set_ylabel(ylabel)
    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(value) for value in light_inputs])
    ax.set_ylim(ylim)
    ax.legend(frameon=False)

    output_path = os.path.join(root_folder, output_filename)
    fig.savefig(output_path, format="svg", dpi=300)
    plt.close(fig)
    return output_path


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Create root-level GFP and mCherry comparison scatter plots, or "
            "regenerate individualized plots when one label is provided."
        )
    )
    parser.add_argument(
        "root_folder",
        help="Root folder containing existing label subfolders.",
    )
    parser.add_argument(
        "-l",
        "--labels",
        nargs="+",
        required=True,
        help=(
            "Label subfolders. One label regenerates that label's plots in its "
            "folder; multiple labels create root-level comparison plots."
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
            "Labels to show in plot titles/legends. Must match --labels count "
            "when provided. Defaults to the folder/strain labels."
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
    outputs = compare_labels(
        root_folder=args.root_folder,
        labels=args.labels,
        num_cond=args.num_cond,
        triplicate=args.triplicate,
        gain=args.gain,
        light_inputs=args.light_inputs,
        gfp_ylim=args.gfp_ylim,
        mcherry_ylim=args.mcherry_ylim,
        plot_labels=args.plot_labels,
    )

    for output_name, output_path in outputs.items():
        if isinstance(output_path, dict):
            print(f"{output_name}:")
            for nested_name, nested_path in output_path.items():
                print(f"  {nested_name}: {nested_path}")
        else:
            print(f"{output_name}: {output_path}")


if __name__ == "__main__":
    main()


__all__ = ["compare_labels"]

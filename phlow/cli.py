import argparse
import os

from .flow_compare import compare_labels
from .flow_pipeline import _parse_light_inputs, _parse_ylim, run_flow_experiment


def _add_common_arguments(parser, labels_required):
    parser.add_argument(
        "--address",
        default=os.getcwd(),
        help="Path to the mother folder. Defaults to the current working directory.",
    )
    parser.add_argument(
        "-l",
        "--labels",
        nargs="+",
        required=labels_required,
        help="Experiment labels / label subfolders.",
    )
    parser.add_argument(
        "-n",
        "--num-cond",
        "--num_cond",
        type=int,
        default=4,
        dest="num_cond",
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
        help="Labels to show in plots/legends. Defaults to folder labels.",
    )
    parser.add_argument(
        "--gfp-ylim",
        type=_parse_ylim,
        help="GFP scatter y-axis limits as min,max. Default: automatic from data.",
    )
    parser.add_argument(
        "--mcherry-ylim",
        type=_parse_ylim,
        help="mCherry scatter y-axis limits as min,max. Default: automatic from data.",
    )


def build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="phlow",
        description="Organize, analyze, and compare flow cytometry .fcs experiments.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run",
        help="Organize raw .fcs files when needed and generate per-label plots.",
    )
    _add_common_arguments(run_parser, labels_required=False)
    run_parser.add_argument(
        "-r",
        "--reverse-numbering",
        action="store_true",
        help="Reverse .fcs numbering inside each label folder after normalization.",
    )
    run_parser.set_defaults(func=_run)

    compare_parser = subparsers.add_parser(
        "compare",
        help="Generate comparison plots for existing label folders.",
    )
    _add_common_arguments(compare_parser, labels_required=True)
    compare_parser.set_defaults(func=_compare)

    return parser


def _print_nested_outputs(outputs):
    for name, output in outputs.items():
        if isinstance(output, dict):
            print(f"{name}:")
            for nested_name, nested_path in output.items():
                print(f"  {nested_name}: {nested_path}")
        else:
            print(f"{name}: {output}")


def _run(args):
    outputs = run_flow_experiment(
        root_folder=args.address,
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
    _print_nested_outputs(outputs)


def _compare(args):
    outputs = compare_labels(
        root_folder=args.address,
        labels=args.labels,
        num_cond=args.num_cond,
        triplicate=args.triplicate,
        gain=args.gain,
        light_inputs=args.light_inputs,
        gfp_ylim=args.gfp_ylim,
        mcherry_ylim=args.mcherry_ylim,
        plot_labels=args.plot_labels,
    )
    _print_nested_outputs(outputs)


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    args.func(args)

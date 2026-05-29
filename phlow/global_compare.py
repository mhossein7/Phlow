import os

from .flow_compare import (
    _labels_filename_suffix,
    _metric_ylim,
    _save_comparison_scatter,
    _save_gfp_histogram_comparison,
)
from .flow_pipeline import (
    DEFAULT_LIGHT_INPUTS,
    _condition_labels,
    _count_fcs_files,
    _label_folder,
    _resolve_plot_labels,
    load_label_flow_unit,
)


GLOBAL_COMPARISON_FOLDER = "Global comparisons"


def _as_list(value, name):
    if value is None:
        return None
    if isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a list")
    return list(value)


def _validate_same_length(values, expected_len, name):
    if len(values) != expected_len:
        raise ValueError(f"{name} must have one value per address")


def _default_plot_labels(addresses, labels):
    plot_labels = []
    for address, label in zip(addresses, labels):
        address_label = os.path.basename(os.path.normpath(address)) or address
        plot_labels.append(f"{label} ({address_label})")
    return plot_labels


def _comparison_folder_name(labels):
    return _labels_filename_suffix(labels).replace("__", "-")


def _global_output_dirs(addresses, labels):
    comparison_folder = _comparison_folder_name(labels)
    output_dirs = []
    for address in addresses:
        output_dir = os.path.join(address, GLOBAL_COMPARISON_FOLDER, comparison_folder)
        os.makedirs(output_dir, exist_ok=True)
        output_dirs.append(output_dir)
    return output_dirs


def global_compare(
    addresses,
    labels,
    num_cond=4,
    gains=None,
    triplicates=None,
    light_inputs=None,
    gfp_ylim=None,
    mcherry_ylim=None,
    plot_labels=None,
    gfp_y_label="Mean log10(GFP)",
    mcherry_y_label="Mean log10(mCherry)",
):
    """
    Compare one label from each experiment folder and save global comparison
    plots into a "Global comparisons" folder under every experiment address.
    """
    addresses = _as_list(addresses, "addresses")
    labels = _as_list(labels, "labels")
    if not addresses or len(addresses) < 2:
        raise ValueError("addresses must contain at least two experiment folders")
    if not labels:
        raise ValueError("labels must contain at least one label")
    _validate_same_length(labels, len(addresses), "labels")

    addresses = [os.path.abspath(str(address)) for address in addresses]
    labels = [str(label) for label in labels]
    for address in addresses:
        if not os.path.isdir(address):
            raise ValueError(f"address does not exist or is not a directory: {address}")

    if gains is None:
        gains = [8] * len(addresses)
    else:
        gains = list(gains)
        _validate_same_length(gains, len(addresses), "gains")

    if triplicates is None:
        triplicates = [False] * len(addresses)
    else:
        triplicates = list(triplicates)
        _validate_same_length(triplicates, len(addresses), "triplicates")

    if light_inputs is None:
        light_inputs = DEFAULT_LIGHT_INPUTS if num_cond == 4 else list(range(1, num_cond + 1))
    else:
        light_inputs = list(light_inputs)

    cond_labels = _condition_labels(num_cond, light_inputs)
    if plot_labels is None:
        plot_labels = _default_plot_labels(addresses, labels)
    comparison_keys = [str(idx) for idx in range(len(addresses))]
    plot_label_map = _resolve_plot_labels(comparison_keys, plot_labels)

    units = {}
    for idx, (address, label, gain, triplicate) in enumerate(
        zip(addresses, labels, gains, triplicates)
    ):
        folder = _label_folder(address, label)
        if not os.path.isdir(folder):
            raise ValueError(f"Label folder does not exist: {folder}")

        files_per_label = num_cond * (3 if triplicate else 1)
        actual_files = _count_fcs_files(folder)
        if actual_files != files_per_label:
            raise ValueError(
                f"Expected {files_per_label} .fcs files in {folder}, found {actual_files}"
            )

        key = comparison_keys[idx]
        units[key] = load_label_flow_unit(
            label=label,
            root_folder=address,
            num_cond=num_cond,
            triplicate=triplicate,
            gain=gain,
            cond_labels=cond_labels,
            plot_label=plot_label_map[key],
        )

    if gfp_ylim is None:
        gfp_ylim = _metric_ylim(units, "GFP")
    if mcherry_ylim is None:
        mcherry_ylim = _metric_ylim(units, "mCherry")

    suffix = _labels_filename_suffix(plot_label_map.values())
    output_dirs = _global_output_dirs(addresses, labels)
    outputs = {}
    for output_dir in output_dirs:
        outputs[output_dir] = {
            "gfp_comparison": _save_comparison_scatter(
                root_folder=output_dir,
                units=units,
                plot_label_map=plot_label_map,
                light_inputs=light_inputs,
                metric_name="GFP",
                ylabel=gfp_y_label,
                ylim=gfp_ylim,
                output_filename=f"global_compare_gfp_vs_light_input_{suffix}.svg",
            ),
            "mcherry_comparison": _save_comparison_scatter(
                root_folder=output_dir,
                units=units,
                plot_label_map=plot_label_map,
                light_inputs=light_inputs,
                metric_name="mCherry",
                ylabel=mcherry_y_label,
                ylim=mcherry_ylim,
                output_filename=f"global_compare_mcherry_vs_light_input_{suffix}.svg",
            ),
            "gfp_histogram_comparison": _save_gfp_histogram_comparison(
                root_folder=output_dir,
                units=units,
                plot_label_map=plot_label_map,
                output_filename=f"global_compare_gfp_histograms_{suffix}.svg",
            ),
        }
    return outputs


__all__ = ["GLOBAL_COMPARISON_FOLDER", "global_compare"]

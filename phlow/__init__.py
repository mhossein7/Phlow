"""Phlow: small tools for flow cytometry .fcs organization and analysis."""

__all__ = ["compare_labels", "flow_unit", "run_flow_experiment"]


def __getattr__(name):
    if name == "compare_labels":
        from .flow_compare import compare_labels

        return compare_labels
    if name == "flow_unit":
        from .flow_object import flow_unit

        return flow_unit
    if name == "run_flow_experiment":
        from .flow_pipeline import run_flow_experiment

        return run_flow_experiment
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialog: create it and re-lay it out for every choice
that drives the layout. Skipped when no display is available."""

import pytest


@pytest.fixture()
def tk_node():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw
    import seamm

    Pmw.initialise(root)
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.tk"
    )
    node = flowchart.create_node("Vibrational Analysis")
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get("Vibrational Analysis")
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    yield tk_node
    root.destroy()


def _mapped(widget):
    return widget.grid_info() != {}


def test_dialog_layouts(tk_node):
    tk_node.create_dialog()
    # Defaults: optimize, with its settings; resonances with their thresholds
    for key in ("hessian route", "optimize", "fmax", "max steps", "derivatives",
                "displacement", "resonances", "resonance window",
                "resonance threshold", "low frequency"):  # fmt: skip
        assert _mapped(tk_node[key]), key

    # No optimization: its settings are hidden
    tk_node["optimize"].set("no")
    tk_node.reset_dialog()
    assert not _mapped(tk_node["fmax"]) and not _mapped(tk_node["max steps"])

    # Separate tasks cannot optimize: the choice itself is hidden and off
    tk_node["optimize"].set("yes")
    tk_node["hessian route"].set("separate tasks")
    tk_node.reset_dialog()
    assert not _mapped(tk_node["optimize"]) and not _mapped(tk_node["fmax"])
    assert tk_node["optimize"].get() == "no"

    # Plain VPT2: no resonance thresholds
    tk_node["resonances"].set("none (plain VPT2)")
    tk_node.reset_dialog()
    assert not _mapped(tk_node["resonance window"])
    assert not _mapped(tk_node["resonance threshold"])
    assert _mapped(tk_node["low frequency"])

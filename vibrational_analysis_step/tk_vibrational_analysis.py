# -*- coding: utf-8 -*-

"""The graphical part of a Vibrational Analysis step"""

import pprint  # noqa: F401
import tkinter as tk  # noqa: F401

import vibrational_analysis_step  # noqa: F401
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


class TkVibrationalAnalysis(seamm.TkNode):
    """
    The graphical part of a Vibrational Analysis step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    namespace : str
        The namespace of the current step.
    tk_subflowchart : TkFlowchart
        A graphical Flowchart representing a subflowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in Vibrational Analysis_parameters.py

    See Also
    --------
    VibrationalAnalysis, TkVibrationalAnalysis,
    VibrationalAnalysisParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )

    def create_dialog(self):
        """Create the dialog: the optimization, the route and derivatives, the
        displacement and the treatment of resonances and low modes."""
        frame = super().create_dialog(title="Vibrational Analysis")
        P = self.node.parameters
        for key in P:
            if key[0] != "_" and key != "results":
                self[key] = P[key].widget(frame)
        # The layout follows these choices
        for key in ("optimize", "hessian route", "resonances"):
            self[key].bind("<<ComboboxSelected>>", self.reset_dialog)
            self[key].bind("<Return>", self.reset_dialog)
            self[key].bind("<FocusOut>", self.reset_dialog)
        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Lay out the widgets, hiding those that do not apply: the optimization
        settings when the structure is not optimized (or the route is separate
        tasks, which cannot optimize), and the resonance thresholds when
        resonances are not treated."""
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        tasks = self["hessian route"].get() == "separate tasks"
        if tasks:
            self["optimize"].set("no")
        optimize = self["optimize"].get() == "yes" and not tasks
        resonances = not self["resonances"].get().startswith("none")

        keys = ["hessian route"]
        if not tasks:
            keys.append("optimize")
        if optimize:
            keys.extend(["fmax", "max steps"])
        keys.extend(["derivatives", "displacement", "resonances"])
        if resonances:
            keys.extend(["resonance window", "resonance threshold"])
        keys.append("low frequency")

        widgets = []
        for row, key in enumerate(keys):
            indent = key in (
                "fmax",
                "max steps",
                "resonance window",
                "resonance threshold",
            )
            self[key].grid(
                row=row,
                column=0,
                columnspan=2 if not indent else 1,
                sticky=tk.EW,
                padx=(30 if indent else 0, 0),
            )
            widgets.append(self[key])
        sw.align_labels(widgets, sticky=tk.E)

        if "results" in self.node.metadata and "results" in self.node.parameters:
            self.setup_results()

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkVibrationalAnalysis.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

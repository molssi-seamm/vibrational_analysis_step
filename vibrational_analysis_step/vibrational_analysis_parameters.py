# -*- coding: utf-8 -*-
"""
Control parameters for the Vibrational Analysis step in a SEAMM flowchart
"""

import logging
import seamm
import pprint  # noqa: F401

logger = logging.getLogger(__name__)


class VibrationalAnalysisParameters(seamm.Parameters):
    """
    The control parameters for Vibrational Analysis.

    You need to replace the "time" entry in dictionary below these comments with the
    definitions of parameters to control this step. The keys are parameters for the
    current plugin,the values are dictionaries as outlined below.

    Examples
    --------
    ::

        parameters = {
            "time": {
                "default": 100.0,
                "kind": "float",
                "default_units": "ps",
                "enumeration": tuple(),
                "format_string": ".1f",
                "description": "Simulation time:",
                "help_text": ("The time to simulate in the dynamics run.")
            },
        }

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
        Each key of the dictionary is a dictionary that contains the
        the following keys:

    parameters["default"] :
        The default value of the parameter, used to reset it.

    parameters["kind"] : enum()
        Specifies the kind of a variable. One of  "integer", "float", "string",
        "boolean", or "enum"

        While the "kind" of a variable might be a numeric value, it may still have
        enumerated custom values meaningful to the user. For instance, if the parameter
        is a convergence criterion for an optimizer, custom values like "normal",
        "precise", etc, might be adequate. In addition, any parameter can be set to a
        variable of expression, indicated by having "$" as the first character in the
        field. For example, $OPTIMIZER_CONV.

    parameters["default_units"] : str
        The default units, used for resetting the value.

    parameters["enumeration"] : tuple
        A tuple of enumerated values.

    parameters["format_string"] : str
        A format string for "pretty" output.

    parameters["description"] : str
        A short string used as a prompt in the GUI.

    parameters["help_text"] : str
        A longer string to display as help for the user.

    See Also
    --------
    VibrationalAnalysis, TkVibrationalAnalysis, VibrationalAnalysisParameters,
    VibrationalAnalysisStep
    """

    parameters = {
        "optimize": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Optimize the structure first:",
            "help_text": (
                "Optimize the structure over the Model Chemistry's engine before the "
                "analysis. VPT2 needs a tight minimum: the cubic and quartic force "
                "constants are differences of Hessians or gradients, and a residual "
                "force shows up in them."
            ),
        },
        "fmax": {
            "default": 0.0005,
            "kind": "float",
            "default_units": "eV/Å",
            "enumeration": tuple(),
            "format_string": ".5f",
            "description": "Converge to forces below:",
            "help_text": "The largest force at which the optimization stops.",
        },
        "max steps": {
            "default": 500,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "d",
            "description": "Maximum optimization steps:",
            "help_text": "The most steps the optimization may take.",
        },
        "hessian route": {
            "default": "automatic",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("automatic", "warm engine over MDI", "separate tasks"),
            "format_string": "",
            "description": "Calculations via:",
            "help_text": (
                "How the Model Chemistry is driven: automatically (the warm MDI engine "
                "when the program has one and is installed here, else separate tasks "
                "where the job runs them); the warm engine, which also allows the "
                "optimization; or separate tasks, which run in parallel but cannot "
                "optimize the structure."
            ),
        },
        "derivatives": {
            "default": "automatic",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("automatic", "Hessians", "gradients"),
            "format_string": "",
            "description": "At the displaced structures:",
            "help_text": (
                "What is computed at the structures displaced along the normal modes: "
                "Hessians (2 per mode) when the engine has an analytic Hessian, else "
                "gradients (about 2n^2 for n modes, but each is cheap and they can run "
                "as parallel tasks). Hessians by finite differences of forces are "
                "possible but cost more than the gradients."
            ),
        },
        "displacement": {
            "default": "automatic",
            "kind": "float",
            "default_units": "",
            "enumeration": ("automatic", "mass-weighted"),
            "format_string": ".3f",
            "description": "Displacement:",
            "help_text": (
                "The step along each normal mode, in dimensionless (reduced) normal "
                "coordinates. 'automatic' is 0.20 for Hessians and 0.10 for "
                "gradients (whose quartic constants use steps up to 3 times as "
                "large); 'mass-weighted' is Gaussian's 0.025 Å amu^1/2 for every "
                "mode. Smaller steps amplify the noise in the Hessians or gradients, "
                "larger ones add higher-order error. Against Gaussian 09 for "
                "propane, 0.05 gave a ZPE 44 cm^-1 off, 0.20 0.2 cm^-1."
            ),
        },
        "resonances": {
            "default": "GVPT2 (variational)",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "GVPT2 (variational)",
                "deperturbed VPT2",
                "none (plain VPT2)",
            ),
            "format_string": "",
            "description": "Fermi resonances:",
            "help_text": (
                "How Fermi resonances are treated: found by Martin's criteria, removed "
                "from the perturbation sums, then the resonant states mixed "
                "variationally (GVPT2, the default); only removed (deperturbed VPT2); "
                "or not treated, which can give wildly wrong frequencies near a "
                "resonance."
            ),
        },
        "resonance window": {
            "default": 200.0,
            "kind": "float",
            "default_units": "1/cm",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Resonance window:",
            "help_text": (
                "Martin's first criterion: how close a fundamental must be to an "
                "overtone or combination band to be checked for a resonance."
            ),
        },
        "resonance threshold": {
            "default": 1.0,
            "kind": "float",
            "default_units": "1/cm",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Resonance threshold:",
            "help_text": (
                "Martin's second criterion: the smallest second-order shift that "
                "counts as a resonance."
            ),
        },
        "low frequency": {
            "default": 200.0,
            "kind": "float",
            "default_units": "1/cm",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Flag modes below:",
            "help_text": (
                "Modes below this frequency are flagged: VPT2 is unreliable for "
                "torsions and other large-amplitude motions. They stay in the "
                "anharmonic zero-point energy."
            ),
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("VibrationalAnalysisParameters.__init__")

        super().__init__(
            defaults={**VibrationalAnalysisParameters.parameters, **defaults}, data=data
        )

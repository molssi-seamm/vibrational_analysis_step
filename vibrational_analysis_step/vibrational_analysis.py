# -*- coding: utf-8 -*-

"""Non-graphical part of the Vibrational Analysis step in a SEAMM flowchart"""

import importlib.resources
import json
import logging
from pathlib import Path
import pprint  # noqa: F401

import numpy as np

import vibrational_analysis_step
from vibrational_analysis_step import vpt2
import molsystem
import seamm
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Vibrational Analysis")

# Add this module's properties to the standard properties
path = importlib.resources.files("vibrational_analysis_step") / "data"
csv_file = path / "properties.csv"
if csv_file.exists():
    molsystem.add_properties_from_file(csv_file)


class VibrationalAnalysis(seamm.Node):
    """
    The non-graphical part of a Vibrational Analysis step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : VibrationalAnalysisParameters
        The control parameters for Vibrational Analysis.

    See Also
    --------
    TkVibrationalAnalysis,
    VibrationalAnalysis, VibrationalAnalysisParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Vibrational Analysis",
        extension=None,
        logger=logger,
    ):
        """A step for Vibrational Analysis in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Vibrational Analysis {self}")

        super().__init__(
            flowchart=flowchart,
            title="Vibrational Analysis",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self._metadata = vibrational_analysis_step.metadata
        self.parameters = vibrational_analysis_step.VibrationalAnalysisParameters()

    @property
    def version(self):
        """The semantic version of this module."""
        return vibrational_analysis_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return vibrational_analysis_step.__git_revision__

    def description_text(self, P=None):
        """A short description of what this step will do."""
        if not P:
            P = self.parameters.values_to_dict()
        text = (
            "Anharmonic vibrational analysis by second-order vibrational perturbation "
            "theory (VPT2) with the flowchart's Model Chemistry."
        )
        if P["optimize"] in (True, "yes"):
            text += f" The structure is first optimized to forces below {P['fmax']}."
        text += (
            f" Then {P['derivatives']} derivatives at structures displaced by "
            f"{P['displacement']} along each normal mode give the cubic and quartic "
            "force constants."
        )
        resonances = str(P["resonances"])
        if resonances.startswith("GVPT2"):
            text += " Fermi resonances are treated variationally (GVPT2)."
        elif resonances.startswith("deperturbed"):
            text += " Fermi resonances are deperturbed."
        else:
            text += " Fermi resonances are not treated."
        text += f" Modes below {P['low frequency']} are flagged."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run the VPT2 analysis of the current structure."""
        next_node = super().run(printer)
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        printer.important(__(self.description_text(P), indent=self.indent))
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        _, configuration = self.get_system_configuration(None)
        result = self.vpt2_analysis(P, configuration, directory)
        data = self.analyze(result=result)

        self.store_results(
            configuration=configuration,
            data=data,
            create_tables=True,
        )
        for key, note in (
            ("mills1972", "The VPT2 formulas."),
            ("barone2005", "Generalized VPT2 (GVPT2) and the Fermi resonances."),
            ("martin1995", "Martin's criteria for Fermi resonances."),
            ("franke2021", "The resonance and polyad treatment."),
            ("pyvpt2", "The VPT2 algebra adapted in this step."),
        ):
            if key in self._bibliography:
                self.references.cite(
                    raw=self._bibliography[key],
                    alias=key,
                    module="vibrational_analysis_step",
                    level=1,
                    note=note,
                )
        return next_node

    def vpt2_analysis(self, P, configuration, directory):
        """Optimize (optionally), take the Hessian, the derivatives at the
        displaced structures and the VPT2 results for the configuration.

        Parameters
        ----------
        P : dict
            The step's parameter values.
        configuration : molsystem._Configuration
            Updated by the optimization.
        directory : pathlib.Path
            Where vpt2.json is written.

        Returns
        -------
        dict
            As :func:`vibrational_analysis_step.vpt2.vpt2`, plus 'hessian',
            'derivatives', 'spread' and 'energy' (kJ/mol).
        """
        from seamm_exec.hessian import HessianDriver

        driver = HessianDriver(
            self,
            route=P["hessian route"],
            who="Vibrational Analysis",
            print_to=printer,
            precision="tight",
        )
        mc = driver.model_chemistry()
        if mc is None:
            raise RuntimeError(
                "Vibrational Analysis needs a Model Chemistry: add a 'Model "
                "Chemistry' step before it to choose the method."
            )
        self.model = mc["level"]
        masses = np.asarray(configuration.atoms.atomic_masses, dtype=float)
        step = float(P["displacement"])
        optimize = P["optimize"] in (True, "yes")
        fmax = P["fmax"].m_as("eV/Å") if hasattr(P["fmax"], "m_as") else P["fmax"]

        with driver.session(configuration) as session:
            if session is None:
                evaluator, _, _, coords = driver.last_plan
                if optimize:
                    printer.important(
                        __(
                            "The structure is not optimized: the Model Chemistry runs "
                            "as separate tasks, with no engine to drive the "
                            "optimization.",
                            indent=self.indent + 4 * " ",
                        )
                    )
                E0, H0, how = driver.fd_hessian_tasks(evaluator, configuration, coords)
                harmonic = self._harmonic(H0, masses, coords)
                derivatives = "gradients"
                phi3, phi4, spread = self._from_gradients(
                    driver, evaluator, configuration, harmonic, step, None
                )
            else:
                if optimize:
                    session.optimize(fmax, int(P["max steps"]))
                E0 = session.energy(units="kJ/mol")
                H0, how = session.hessian()
                coords = session.coordinates
                harmonic = self._harmonic(H0, masses, coords)
                derivatives = P["derivatives"]
                if derivatives == "automatic":
                    derivatives = "Hessians" if session.analytic else "gradients"
                if derivatives == "Hessians":
                    if not session.analytic:
                        printer.important(
                            __(
                                "The engine has no analytic Hessian, so each displaced "
                                "Hessian is a finite difference of its forces.",
                                indent=self.indent + 4 * " ",
                            )
                        )
                    # The reference Hessian is the one above, from the same
                    # engine and settings as the displaced ones
                    hessians = {
                        key: session.hessian(xyz * vpt2.BOHR_TO_ANGSTROM)[0]
                        for key, xyz in vpt2.displacements(harmonic, "hessians", step)
                    }
                    phi3, phi4, spread = vpt2.force_constants_from_hessians(
                        harmonic, hessians, step, reference=H0
                    )
                else:
                    phi3, phi4, spread = self._from_gradients(
                        driver,
                        session.evaluator,
                        configuration,
                        harmonic,
                        step,
                        session,
                    )

        resonances = str(P["resonances"])
        result = vpt2.vpt2(
            harmonic,
            phi3,
            phi4,
            fermi=not resonances.startswith("none"),
            gvpt2=resonances.startswith("GVPT2"),
            omega_threshold=self._in(P["resonance window"], "1/cm"),
            k_threshold=self._in(P["resonance threshold"], "1/cm"),
            low_mode=self._in(P["low frequency"], "1/cm"),
        )
        result["hessian"] = how
        result["derivatives"] = derivatives
        result["spread"] = spread
        result["energy"] = E0
        self._write_files(directory, harmonic, phi3, phi4, result)
        return result

    #: The spreads between equivalent finite-difference routes (cm^-1) above
    #: which the force constants are too noisy to trust to a few cm^-1
    NOISE = {"cubic": 2.0, "quartic": 5.0}

    @classmethod
    def noisy(cls, spread):
        """Whether the force constants are too noisy (see ``NOISE``)."""
        return any(spread.get(key, 0.0) > limit for key, limit in cls.NOISE.items())

    @staticmethod
    def _in(value, units):
        """A parameter's magnitude in the given units."""
        return float(value.m_as(units)) if hasattr(value, "m_as") else float(value)

    def _harmonic(self, H0, masses, coords_ang):
        """The harmonic analysis; a clear error for molecules VPT2 here cannot
        treat."""
        try:
            return vpt2.Harmonic(
                H0, masses, np.asarray(coords_ang) / vpt2.BOHR_TO_ANGSTROM
            )
        except (vpt2.DegenerateModesError, vpt2.NotAMinimumError) as e:
            raise RuntimeError(f"Vibrational Analysis: {e}") from None

    def _from_gradients(
        self, driver, evaluator, configuration, harmonic, step, session
    ):
        """The force constants from gradients at the displaced structures: as
        tasks (in parallel) when the evaluator runs the model chemistry as
        tasks, otherwise over the warm engine."""
        geometries = {
            "ref": harmonic.coordinates * vpt2.BOHR_TO_ANGSTROM,
        }
        keys = {}
        for key, xyz in vpt2.displacements(harmonic, "gradients", step):
            name = "d" + "_".join(f"{m}{'p' if s > 0 else 'm'}{abs(s)}" for m, s in key)
            keys[name] = key
            geometries[name] = xyz * vpt2.BOHR_TO_ANGSTROM
        printer.important(
            __(
                f"{len(geometries)} gradients for the force constants.",
                indent=self.indent + 4 * " ",
            )
        )
        if session is None or (
            evaluator.path == "batch" and not driver.route.startswith("warm engine")
        ):
            results = driver.evaluate(
                evaluator, configuration, geometries, "the force constants"
            )
            gradients = {name: g for name, (E, g) in results.items()}
        else:
            gradients = {name: session.gradient(x) for name, x in geometries.items()}
        g0 = gradients.pop("ref")
        by_key = {keys[name]: g for name, g in gradients.items()}
        return vpt2.force_constants_from_gradients(harmonic, by_key, step, reference=g0)

    def _write_files(self, directory, harmonic, phi3, phi4, result):
        """The force constants and results, for reuse and checking."""
        data = {
            "masses": harmonic.masses.tolist(),
            "coordinates (bohr)": harmonic.coordinates.tolist(),
            "omega": harmonic.omega.tolist(),
            "rotational constants": harmonic.B.tolist(),
            "phi_ijk": phi3.tolist(),
            "phi_iijj": phi4.tolist(),
        }
        for key in ("nu", "nu_deperturbed", "overtones", "combinations", "chi"):
            data[key] = np.asarray(result[key]).tolist()
        for key in ("zpe_harmonic", "zpe", "g0", "hessian", "derivatives", "spread"):
            data[key] = result[key]
        data["resonances"] = [
            {
                k: (float(v) if isinstance(v, (np.floating, float)) else v)
                for k, v in r.items()
            }
            for r in result["resonances"]
        ]
        data["flagged"] = result["flagged"]
        (directory / "vpt2.json").write_text(json.dumps(data, indent=2, default=str))

    def analyze(self, indent="", result=None, **kwargs):
        """Print the results and return them for storing."""
        if result is None:
            return {}
        cm_to_kj = vpt2.CM_TO_KJMOL
        omega = np.asarray(result["omega"])
        nu = np.asarray(result["nu"])
        flagged = set(result["flagged"])
        resonant = {r["fundamental"] for r in result["resonances"]}
        lines = [
            "",
            "    Mode   Harmonic   Anharmonic      Shift   Overtone",
            "             (cm^-1)      (cm^-1)    (cm^-1)    (cm^-1)",
            "    ----   --------   ----------   --------   --------",
        ]
        for i, (w, v, o) in enumerate(zip(omega, nu, result["overtones"])):
            notes = []
            if i in resonant:
                notes.append("Fermi resonance")
            if i in flagged:
                notes.append("low frequency: VPT2 unreliable")
            lines.append(
                f"    {i + 1:4d}   {w:8.2f}   {v:10.2f}   {v - w:8.2f}   {o:8.2f}"
                + (f"   ({'; '.join(notes)})" if notes else "")
            )
        zpe_h, zpe = result["zpe_harmonic"], result["zpe"]
        text = (
            f"The harmonic Hessian ({result['hessian']}) and the "
            f"{result['derivatives']} at the displaced structures give these "
            "anharmonic frequencies:"
        )
        printer.normal(__(text, indent=self.indent + 4 * " "))
        printer.normal("\n".join(lines))
        printer.normal("")
        text = (
            f"Zero-point energy: harmonic {zpe_h * cm_to_kj:.3f} kJ/mol, anharmonic "
            f"{zpe * cm_to_kj:.3f} kJ/mol (correction {(zpe - zpe_h) * cm_to_kj:.3f} "
            f"kJ/mol, {zpe - zpe_h:.2f} cm^-1, including G0 = {result['g0']:.2f} "
            "cm^-1)."
        )
        printer.normal(__(text, indent=self.indent + 4 * " "))
        if result["resonances"]:
            items = []
            for r in result["resonances"]:
                i, j = r["states"]
                k = r["fundamental"]
                band = f"2 nu{i + 1}" if r["type"] == 1 else f"nu{i + 1} + nu{j + 1}"
                items.append(
                    f"nu{k + 1} ~ {band} (mismatch {r['delta']:.1f} cm^-1, K = "
                    f"{r['K']:.1f} cm^-1)"
                )
            printer.normal(
                __(
                    "Fermi resonances treated: " + "; ".join(items) + ".",
                    indent=self.indent + 4 * " ",
                )
            )
        if flagged:
            printer.normal(
                __(
                    "Modes "
                    + ", ".join(str(i + 1) for i in sorted(flagged))
                    + " are low-frequency (torsions or other large-amplitude "
                    "motions): VPT2 is unreliable for them, though they are kept in "
                    "the zero-point energy.",
                    indent=self.indent + 4 * " ",
                )
            )
        spread = result["spread"]
        printer.normal(
            __(
                "Numerical noise: equivalent finite-difference routes to the force "
                f"constants differ by up to {spread['cubic']:.2f} cm^-1 (cubic) and "
                f"{spread['quartic']:.2f} cm^-1 (quartic).",
                indent=self.indent + 4 * " ",
            )
        )
        if self.noisy(spread):
            printer.important(
                __(
                    "WARNING: the force constants are noisy, so the anharmonic "
                    "results above may be off by several cm^-1 or more. The "
                    "Hessians or gradients are not precise enough for the finite "
                    "differences: tighten the method's convergence, use analytic "
                    "Hessians if the method has them, or try a larger displacement.",
                    indent=self.indent + 4 * " ",
                )
            )
        printer.normal("")
        return {
            "energy": result["energy"],
            "ZPE": zpe * cm_to_kj,
            "harmonic ZPE": zpe_h * cm_to_kj,
            "ZPE correction": (zpe - zpe_h) * cm_to_kj,
            "G0": result["g0"],
            "harmonic frequencies": omega.tolist(),
            "fundamentals": nu.tolist(),
            "overtones": np.asarray(result["overtones"]).tolist(),
            "combination bands": np.asarray(result["combinations"]).tolist(),
            "chi": np.asarray(result["chi"]).tolist(),
            "resonances": "; ".join(
                f"{r['fundamental'] + 1}:{r['states'][0] + 1},{r['states'][1] + 1}"
                for r in result["resonances"]
            ),
            "flagged modes": [i + 1 for i in sorted(flagged)],
            "n_modes": len(omega),
        }

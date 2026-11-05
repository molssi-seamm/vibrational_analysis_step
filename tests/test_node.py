#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The step's analysis over a fake MDI engine on the model water of
test_vpt2: optimization, the Hessian and the two routes for the displaced
structures, through seamm_exec's shared HessianDriver."""

import contextlib
import json
import types

import numpy as np
import pytest

from seamm_exec.hessian import EngineSession, HessianDriver
from vibrational_analysis_step import VibrationalAnalysis, vpt2

from . import test_vpt2 as model

BOHR = vpt2.BOHR_TO_ANGSTROM
EV = 27.211386245988  # eV per hartree


class ModelEngine:
    """An MDI engine for the model water: energies, forces and (optionally)
    an analytic Hessian."""

    def __init__(self, analytic):
        self._analytic = analytic
        self.x = None  # bohr
        self.n_hessians = 0

    def set_coordinates(self, x, units="Å"):
        x = np.asarray(x, dtype=float).reshape(-1)
        self.x = x / BOHR if units == "Å" else x

    def _energy(self):
        x = self.x.reshape(3, 3)
        r1 = np.linalg.norm(x[1] - x[0])
        r2 = np.linalg.norm(x[2] - x[0])
        e = sum(
            model.D * (1 - np.exp(-model.A * (r - model.RE))) ** 2 for r in (r1, r2)
        )
        e += model.C_RR * (r1 - model.RE) * (r2 - model.RE)
        cos = (x[1] - x[0]) @ (x[2] - x[0]) / (r1 * r2)
        e += 0.5 * model.K_BEND * (np.arccos(np.clip(cos, -1, 1)) - model.THETA_E) ** 2
        return e

    def energy(self, units="kJ/mol"):
        e = self._energy()
        return {"kJ/mol": e * 2625.4996394799, "eV": e * EV}.get(units, e)

    def forces(self, units="hartree/bohr"):
        f = -model.model_gradient(self.x)
        if units == "eV/Å":
            return (f * EV / BOHR).reshape(-1, 3)
        return f

    def supports(self, command):
        return self._analytic and command == "<HESSIAN"

    def hessian(self):
        self.n_hessians += 1
        return model.model_hessian(self.x)


@pytest.fixture()
def configuration():
    from molsystem import SystemDB

    db = SystemDB(filename="file:vib_test?mode=memory&cache=shared")
    system = db.create_system(name="water")
    c = system.create_configuration(name="model")
    xyz = model.model_minimum() * BOHR
    # Start a little away from the minimum so the optimization has work to do
    xyz[1] += [0.01, -0.005, 0.0]
    c.atoms.append(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], symbol=["O", "H", "H"])
    c.charge = 0
    c.spin_multiplicity = 1
    yield c
    db.close()


def run_analysis(configuration, tmp_path, monkeypatch, analytic, **values):
    node = VibrationalAnalysis()
    node._id = ("1",)
    mc = {"level": "Model:water", "step": "Model", "method": "water", "options": {}}
    monkeypatch.setattr(
        node, "variable_exists", lambda name: name == "_model_chemistry"
    )
    monkeypatch.setattr(node, "get_variable", lambda name: mc)
    engine = ModelEngine(analytic)

    @contextlib.contextmanager
    def session(self, configuration, announce=True):
        xyz = configuration.atoms.get_coordinates(fractionals=False, as_array=True)
        evaluator = types.SimpleNamespace(path="mdi")
        yield EngineSession(
            self, engine, configuration, evaluator, [8, 1, 1], np.asarray(xyz)
        )

    monkeypatch.setattr(HessianDriver, "session", session)
    P = node.parameters.current_values_to_dict(context={})
    P.update(values)
    # The masses of the model test (the most abundant isotopes)
    monkeypatch.setattr(
        type(configuration.atoms),
        "atomic_masses",
        property(lambda self: list(model.MASSES)),
    )
    result = node.vpt2_analysis(P, configuration, tmp_path)
    return node, engine, result


def test_hessian_route(configuration, tmp_path, monkeypatch):
    """Analytic Hessians: the optimized structure, the Hessian route, 6
    displaced Hessians, and the same numbers as the core on the model."""
    node, engine, result = run_analysis(configuration, tmp_path, monkeypatch, True)
    assert result["derivatives"] == "Hessians" and result["hessian"] == "analytic"
    assert engine.n_hessians == 1 + 2 * 3
    # The core directly, at the exact minimum
    xyz = model.model_minimum()
    harmonic = vpt2.Harmonic(model.model_hessian(xyz), model.MASSES, xyz)
    data = {k: model.model_hessian(x) for k, x in vpt2.displacements(harmonic)}
    phi3, phi4, _ = vpt2.force_constants_from_hessians(harmonic, data)
    ref = vpt2.vpt2(harmonic, phi3, phi4)
    assert np.allclose(result["omega"], ref["omega"], atol=0.05)
    assert np.allclose(result["nu"], ref["nu"], atol=0.1)
    assert result["zpe"] == pytest.approx(ref["zpe"], abs=0.05)
    # The optimization moved the structure back to the minimum
    xyz_now = configuration.atoms.get_coordinates(fractionals=False, as_array=True)
    assert np.abs(model.model_gradient(np.asarray(xyz_now) / BOHR)).max() < 1e-4
    saved = json.loads((tmp_path / "vpt2.json").read_text())
    assert saved["derivatives"] == "Hessians" and len(saved["phi_ijk"]) == 3


def test_gradient_route_agrees(configuration, tmp_path, monkeypatch):
    """No analytic Hessian: finite-difference Hessian over the engine, and
    gradients at the displaced structures; the same results."""
    _, engine, hes = run_analysis(configuration, tmp_path, monkeypatch, True)
    _, engine, grad = run_analysis(configuration, tmp_path, monkeypatch, False)
    assert grad["derivatives"] == "gradients"
    assert grad["hessian"] == "finite differences over MDI" and engine.n_hessians == 0
    assert np.allclose(grad["nu"], hes["nu"], atol=0.5)
    assert grad["zpe"] == pytest.approx(hes["zpe"], abs=0.3)


def test_analyze_and_results(configuration, tmp_path, monkeypatch):
    node, _, result = run_analysis(configuration, tmp_path, monkeypatch, True)
    data = node.analyze(result=result)
    assert data["n_modes"] == 3 and len(data["fundamentals"]) == 3
    assert data["ZPE"] < data["harmonic ZPE"]
    assert data["ZPE correction"] == pytest.approx(data["ZPE"] - data["harmonic ZPE"])
    assert data["flagged modes"] == []
    # Every result is described in the metadata
    assert set(data) <= set(node.metadata["results"])


def test_degenerate_molecule_is_refused(tmp_path):
    """A spherical top stops with a clear message (Paul, 2026-11-05)."""
    node = VibrationalAnalysis()
    r = 1.09 / np.sqrt(3)
    xyz = [[0, 0, 0], [r, r, r], [-r, -r, r], [-r, r, -r], [r, -r, -r]]
    with pytest.raises(RuntimeError, match="spherical top"):
        node._harmonic(np.eye(15), [12.0, 1.008, 1.008, 1.008, 1.008], xyz)

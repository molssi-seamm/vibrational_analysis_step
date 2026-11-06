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
    step = node.STEP["Hessians"]
    data = {
        k: model.model_hessian(x) for k, x in vpt2.displacements(harmonic, step=step)
    }
    phi3, phi4, _ = vpt2.force_constants_from_hessians(harmonic, data, step)
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
    # Different default steps (0.20 and 0.10): the higher-order (truncation)
    # error of this strongly anharmonic model differs by up to ~0.75 cm^-1
    assert np.allclose(grad["nu"], hes["nu"], atol=1.0)
    assert grad["zpe"] == pytest.approx(hes["zpe"], abs=0.5)


def test_analyze_and_results(configuration, tmp_path, monkeypatch):
    node, _, result = run_analysis(configuration, tmp_path, monkeypatch, True)
    data = node.analyze(result=result)
    assert data["n_modes"] == 3 and len(data["fundamentals"]) == 3
    assert data["ZPE"] < data["harmonic ZPE"]
    assert data["ZPE correction"] == pytest.approx(data["ZPE"] - data["harmonic ZPE"])
    assert data["flagged modes"] == []
    # Every result is described in the metadata
    assert set(data) <= set(node.metadata["results"])


def _tetrahedral_springs():
    """A tetrahedral CH4-like molecule held by springs between every pair of
    atoms: a spherical top with doubly and triply degenerate modes."""
    r = 1.09 / np.sqrt(3) / BOHR
    xyz = np.array([[0, 0, 0], [r, r, r], [-r, -r, r], [-r, r, -r], [r, -r, -r]])
    n = len(xyz)
    H = np.zeros((3 * n, 3 * n))
    for a in range(n):
        for b in range(a + 1, n):
            d = xyz[b] - xyz[a]
            u = d / np.linalg.norm(d)
            block = 0.3 * np.outer(u, u)
            for i, j, sign in ((a, a, 1), (b, b, 1), (a, b, -1), (b, a, -1)):
                H[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] += sign * block
    return H, xyz


def test_spherical_top_gives_zpe_only():
    """A spherical top is accepted for the ZPE only (Paul, Stage A, 2026-11-05):
    no anharmonic fundamentals are reported."""
    node = VibrationalAnalysis()
    node._id = ("1",)
    H, xyz = _tetrahedral_springs()
    masses = [12.0, 1.00782503, 1.00782503, 1.00782503, 1.00782503]
    harmonic = node._harmonic(H, masses, xyz * BOHR)
    assert harmonic.rotor == "spherical top" and harmonic.degenerate
    assert [len(s) for s in harmonic.degenerate_sets] == [3, 2, 3]
    n = harmonic.n
    result = vpt2.vpt2(harmonic, np.zeros((n, n, n)), np.zeros((n, n)))
    assert result["zpe_only"]
    result.update(
        hessian="analytic", derivatives="Hessians", energy=0.0,
        spread={"cubic": 0.0, "quartic": 0.0},
    )  # fmt: skip
    data = node.analyze(result=result)
    assert "fundamentals" not in data and "overtones" not in data
    assert data["ZPE"] == pytest.approx(data["harmonic ZPE"] + data["ZPE correction"])
    assert set(data) <= set(node.metadata["results"])


def test_linear_molecule_is_refused():
    """Linear molecules (3N-5 modes) are still refused."""
    node = VibrationalAnalysis()
    xyz = [[0, 0, -1.16], [0, 0, 0], [0, 0, 1.16]]
    with pytest.raises(RuntimeError, match="linear"):
        node._harmonic(np.eye(9), [15.995, 12.0, 15.995], xyz)


def test_noise_gate():
    assert not VibrationalAnalysis.noisy({"cubic": 0.3, "quartic": 1.9})
    # A large cubic spread at the default steps is truncation, not noise
    assert not VibrationalAnalysis.noisy({"cubic": 8.0, "quartic": 1.0})
    assert VibrationalAnalysis.noisy({"cubic": 0.2, "quartic": 92.0})


def test_default_steps():
    node = VibrationalAnalysis()
    assert node._step("automatic", None, "Hessians") == 0.20
    assert node._step("automatic", None, "gradients") == 0.10
    assert node._step(0.07, None, "Hessians") == 0.07
    h = types.SimpleNamespace(omega=np.array([200.0, 3000.0]))
    steps = node._step("mass-weighted", h, "Hessians")
    assert steps[0] < steps[1] and steps[1] == pytest.approx(0.236, abs=0.005)

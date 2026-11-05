#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of the VPT2 core (vibrational_analysis_step.vpt2 and .fermi).

The Fermi-solver references are pyVPT2's (BSD-3-Clause; see
THIRD_PARTY_LICENSES.rst).
"""

from types import SimpleNamespace

import numpy as np
import pytest

from vibrational_analysis_step import vpt2
from vibrational_analysis_step.fermi import Interaction, Polyad, State, fermi_solver

ANGSTROM = 1 / 0.529177210903  # bohr


# --------------------------------------------------------------------------
# A model water: Morse O-H bonds, a harmonic bend and a stretch-stretch term,
# with an analytic gradient and a Hessian by central differences of it.
# --------------------------------------------------------------------------
D, A, RE = 0.20, 1.2, 1.81
K_BEND, THETA_E = 0.16, np.radians(104.5)
C_RR = 0.01
MASSES = np.array([15.99491462, 1.00782503, 1.00782503])


def model_gradient(xyz):
    xyz = np.asarray(xyz, dtype=float).reshape(3, 3)
    g = np.zeros((3, 3))
    d1, d2 = xyz[1] - xyz[0], xyz[2] - xyz[0]
    r1, r2 = np.linalg.norm(d1), np.linalg.norm(d2)
    u1, u2 = d1 / r1, d2 / r2
    for k, (r, u) in enumerate(((r1, u1), (r2, u2)), start=1):
        e = np.exp(-A * (r - RE))
        dV = 2 * D * (1 - e) * A * e
        g[k] += dV * u
        g[0] -= dV * u
    # Stretch-stretch coupling
    g[1] += C_RR * (r2 - RE) * u1
    g[2] += C_RR * (r1 - RE) * u2
    g[0] -= C_RR * ((r2 - RE) * u1 + (r1 - RE) * u2)
    # Bend
    cos = np.clip(u1 @ u2, -1, 1)
    theta = np.arccos(cos)
    sin = np.sin(theta)
    dV = K_BEND * (theta - THETA_E)
    dt1 = -(u2 - cos * u1) / (r1 * sin)
    dt2 = -(u1 - cos * u2) / (r2 * sin)
    g[1] += dV * dt1
    g[2] += dV * dt2
    g[0] -= dV * (dt1 + dt2)
    return g.ravel()


def model_hessian(xyz, h=1.0e-5):
    x = np.asarray(xyz, dtype=float).ravel()
    H = np.zeros((9, 9))
    for i in range(9):
        p, m = x.copy(), x.copy()
        p[i] += h
        m[i] -= h
        H[:, i] = (model_gradient(p) - model_gradient(m)) / (2 * h)
    return 0.5 * (H + H.T)


def model_minimum():
    """The minimum of the model, in a general orientation, by Newton steps
    (the Hessian's pseudo-inverse skips translations and rotations)."""
    x = np.array([[0.1, -0.05, 0.2], [1.85, 0.1, 0.0], [-0.45, 1.75, 0.15]]).ravel()
    for _ in range(50):
        g = model_gradient(x)
        if np.abs(g).max() < 1e-12:
            break
        x = x - np.linalg.pinv(model_hessian(x), rcond=1e-8) @ g
    assert np.abs(model_gradient(x)).max() < 1e-10
    return x.reshape(3, 3)


@pytest.fixture(scope="module")
def water():
    xyz = model_minimum()
    harmonic = vpt2.Harmonic(model_hessian(xyz), MASSES, xyz)
    return xyz, harmonic


def test_harmonic_water(water):
    xyz, harmonic = water
    assert harmonic.n == 3
    assert harmonic.rotor == "asymmetric top"
    assert np.all(np.diff(harmonic.omega) > 0)
    # Unit step in q_i: the Hessian in normal coordinates is diag(omega)
    Hq = harmonic.to_normal(model_hessian(xyz))
    assert np.allclose(Hq, np.diag(harmonic.omega), atol=1e-3)
    assert harmonic.zpe == pytest.approx(0.5 * harmonic.omega.sum())


def test_displacement_counts(water):
    _, harmonic = water
    assert len(vpt2.displacements(harmonic, "hessians")) == 2 * 3
    assert len(vpt2.displacements(harmonic, "gradients")) == 2 * 9 + 2 * 3


def run_route(harmonic, route, step=0.05):
    geometries = vpt2.displacements(harmonic, route, step)
    if route == "hessians":
        data = {key: model_hessian(x) for key, x in geometries}
        phi3, phi4, spread = vpt2.force_constants_from_hessians(harmonic, data, step)
    else:
        data = {key: model_gradient(x) for key, x in geometries}
        phi3, phi4, spread = vpt2.force_constants_from_gradients(harmonic, data, step)
    return phi3, phi4, spread, vpt2.vpt2(harmonic, phi3, phi4)


def test_routes_agree(water):
    """Hessians and gradients give the same force constants and results."""
    _, harmonic = water
    h3, h4, hspread, hres = run_route(harmonic, "hessians")
    g3, g4, gspread, gres = run_route(harmonic, "gradients")
    assert np.max(np.abs(h3 - g3)) < 1.0
    assert np.max(np.abs(np.diag(h4) - np.diag(g4))) < 1.0
    assert np.allclose(hres["nu"], gres["nu"], atol=0.5)
    assert hres["zpe"] == pytest.approx(gres["zpe"], abs=0.2)
    # The finite-difference routes within each agree with each other
    assert hspread["cubic"] < 1.0 and hspread["quartic"] < 1.0


def test_water_is_anharmonic(water):
    _, harmonic = water
    _, _, _, result = run_route(harmonic, "hessians")
    shift = result["nu"] - result["omega"]
    assert np.all(shift < 0)
    # The stretches drop most
    assert shift[2] < shift[0] and shift[1] < shift[0]
    assert result["zpe"] < result["zpe_harmonic"]
    assert result["flagged"] == []


def test_morse_is_exact():
    """Uncoupled Morse oscillators: VPT2 gives the exact levels.

    In dimensionless coordinates phi_3 = -6 D x^1.5, phi_4 = 14 D x^2 with
    omega = 2 D x; the exact fundamental is omega - 2 omega x_e and the exact
    zero-point energy omega / 2 - omega x_e / 4, with omega x_e = omega^2 / 4D.
    """
    omega = np.array([1100.0, 1700.0, 3900.0])
    Dm = np.array([30000.0, 45000.0, 38000.0])
    x = omega / (2 * Dm)
    n = 3
    phi3 = np.zeros((n, n, n))
    phi4 = np.zeros((n, n))
    for i in range(n):
        phi3[i, i, i] = -6 * Dm[i] * x[i] ** 1.5
        phi4[i, i] = 14 * Dm[i] * x[i] ** 2
    harmonic = SimpleNamespace(
        omega=omega, B=np.zeros(3), zeta=np.zeros((3, n, n)), zpe=omega.sum() / 2
    )
    result = vpt2.vpt2(harmonic, phi3, phi4)
    wx = omega**2 / (4 * Dm)
    assert np.allclose(result["nu"], omega - 2 * wx, atol=1e-9)
    assert np.allclose(np.diag(result["chi"]), -wx, atol=1e-9)
    assert result["zpe"] == pytest.approx(np.sum(omega / 2 - wx / 4), abs=1e-9)
    assert result["g0"] == pytest.approx(0.0, abs=1e-9)
    assert np.allclose(result["overtones"], 2 * omega - 6 * wx, atol=1e-9)


def test_symmetric_top_refused():
    """CH4 (a spherical top) is refused: degenerate modes are not handled."""
    r = 1.09 * ANGSTROM / np.sqrt(3)
    xyz = np.array([[0, 0, 0], [r, r, r], [-r, -r, r], [-r, r, -r], [r, -r, -r]])
    masses = [12.0, 1.00782503, 1.00782503, 1.00782503, 1.00782503]
    with pytest.raises(vpt2.DegenerateModesError, match="spherical top"):
        vpt2.Harmonic(np.eye(15), masses, xyz)


def test_not_a_minimum(water):
    xyz, _ = water
    H = model_hessian(xyz)
    with pytest.raises(vpt2.NotAMinimumError):
        vpt2.Harmonic(-H, MASSES, xyz)


# --------------------------------------------------------------------------
# Fermi resonances; references from pyVPT2's tests
# --------------------------------------------------------------------------
def test_fermi_solver():
    ref = {(1, 2): 3047.2, (0,): 3116.3, (3, 6): 2711.5, (5,): 2851.5, (2, 6): 3004.3}
    interactions = [
        Interaction(State((0,), 3095.5), State((1, 2), 3068.0), -89.7, 2),
        Interaction(State((5,), 2828.0), State((2, 6), 2987.4), -146.7, 2),
        Interaction(State((5,), 2828.0), State((3, 6), 2751.9), 185.8, 2),
    ]
    result = fermi_solver(interactions)
    assert set(result) == set(ref)
    for state, value in ref.items():
        assert result[state] == pytest.approx(value, abs=0.1)


def test_fermi_solver_empty():
    assert fermi_solver([]) == {}


def test_polyad_type1():
    ref = {(0,): 3097.2, (1, 1): 3160.8}
    polyad = Polyad(Interaction(State((0,), 3106.1), State((1, 1), 3152.0), 88.1, 1))
    for state, value in polyad.solve().items():
        assert value == pytest.approx(ref[state], abs=0.1)


def test_find_resonances_type1():
    omega = [1544.3, 1601.8, 2956.6]
    phi = np.zeros((3, 3, 3))
    for (i, j, k), v in {
        (0, 0, 0): 158.9,
        (0, 0, 1): -247.3,
        (0, 0, 2): -166.7,
        (0, 1, 1): 120.1,
        (0, 1, 2): -142.5,
        (0, 2, 2): -281.9,
        (1, 1, 1): -426.8,
        (1, 1, 2): 52.4,
        (1, 2, 2): 151.6,
        (2, 2, 2): 2406.3,
    }.items():
        for p in set(__import__("itertools").permutations((i, j, k))):
            phi[p] = v
    found = vpt2.find_resonances(omega, phi)
    assert [(r["type"], r["fundamental"], r["states"]) for r in found] == [
        (1, 2, (0, 0))
    ]


# --------------------------------------------------------------------------
# Against ORCA 6.1.1's own VPT2: H2O B3LYP/def2-TZVP AutoAux DEFGRID3 TightSCF
# NoCOSX, Hessians from ORCA at the core's displacements (ChemAI 2026-11-05)
# --------------------------------------------------------------------------
def test_h2o_against_orca():
    from pathlib import Path

    data = np.load(Path(__file__).parent / "data" / "h2o_b3lyp_def2tzvp_orca.npz")
    harmonic = vpt2.Harmonic(data["hessian0"], data["masses"], data["coordinates"])
    # ORCA's harmonic frequencies
    assert np.allclose(harmonic.omega, [1616.357, 3782.812, 3888.050], atol=0.02)
    # The fixture's Hessians by their geometries, matched to this machine's
    # displacements (the phases of the modes there may differ)
    hessians = {}
    for key, xyz in vpt2.displacements(harmonic, "hessians", float(data["step"])):
        distance = np.abs(data["geometries"] - xyz[None]).max(axis=(1, 2))
        assert distance.min() < 1e-8
        hessians[key] = data["displaced"][np.argmin(distance)]
    # The reference Hessian must come from the same runs (frame, grid) as the
    # displaced ones: the quartic constants are second differences over 0.05^2
    phi3, phi4, spread = vpt2.force_constants_from_hessians(
        harmonic, hessians, float(data["step"]), reference=data["reference"]
    )
    assert spread["quartic"] < 3.0
    result = vpt2.vpt2(harmonic, phi3, phi4)
    assert np.allclose(result["nu"], [1566.641, 3615.203, 3706.130], atol=0.05)
    # ORCA's "Anharmonic correction" to the ZPE leaves out the constant term
    correction = result["zpe"] - result["g0"] - result["zpe_harmonic"]
    assert correction == pytest.approx(-73.224, abs=0.05)
    assert result["resonances"] == []


def test_type1_resonance_is_deperturbed_in_off_diagonal_chi():
    """omega_2 ~ 2 omega_1: the resonant denominator enters chi_12 through its
    k = 1 term and must be removed there too (propane against Gaussian,
    2026-11-05: chi was +2913 instead of -9.8 cm^-1 at a 0.26 cm^-1 mismatch).
    The deperturbed chi_12 stays finite and smooth through the resonance."""
    n = 2

    def chi12(w2, fermi):
        omega = np.array([1500.0, w2])
        phi3 = np.zeros((n, n, n))
        for p in {(0, 0, 1), (0, 1, 0), (1, 0, 0)}:
            phi3[p] = 80.0
        phi4 = np.diag([30.0, 60.0])
        harmonic = SimpleNamespace(
            omega=omega, B=np.zeros(3), zeta=np.zeros((3, n, n)), zpe=omega.sum() / 2
        )
        r = vpt2.vpt2(harmonic, phi3, phi4, fermi=fermi, gvpt2=False)
        return r["chi"][0, 1], r["zpe"]

    a, za = chi12(3000.1, True)
    b, zb = chi12(3000.3, True)
    assert abs(a) < 50 and abs(a - b) < 0.5 and abs(za - zb) < 0.5
    plain, _ = chi12(3000.1, False)
    assert abs(plain) > 1000


def test_per_mode_steps(water):
    """Per-mode steps (Gaussian's mass-weighted convention) give the same
    force constants as a uniform step, by both routes."""
    _, harmonic = water
    steps = vpt2.mass_weighted_steps(harmonic)
    # Larger dimensionless steps for stiffer modes; 0.025 A amu^1/2 ~ 0.24 at
    # 3000 cm^-1
    assert np.all(np.diff(steps) > 0)
    assert vpt2.mass_weighted_steps(
        SimpleNamespace(omega=np.array([3000.0]))
    ) == pytest.approx(0.236, abs=0.005)
    ref = run_route(harmonic, "hessians", 0.05)[3]
    for route in ("hessians", "gradients"):
        geometries = vpt2.displacements(harmonic, route, steps)
        if route == "hessians":
            data = {k: model_hessian(x) for k, x in geometries}
            f3, f4, _ = vpt2.force_constants_from_hessians(harmonic, data, steps)
        else:
            data = {k: model_gradient(x) for k, x in geometries}
            f3, f4, _ = vpt2.force_constants_from_gradients(harmonic, data, steps)
        r = vpt2.vpt2(harmonic, f3, f4)
        assert np.allclose(r["nu"], ref["nu"], atol=1.0), route
        assert r["zpe"] == pytest.approx(ref["zpe"], abs=0.5)

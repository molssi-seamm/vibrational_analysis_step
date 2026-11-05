# -*- coding: utf-8 -*-

"""Second-order vibrational perturbation theory (VPT2), independent of any code.

The anharmonic frequencies and zero-point energy of a molecule at a minimum, from

1. its Cartesian Hessian there (:class:`Harmonic`: the normal modes and harmonic
   frequencies, rotational constants and Coriolis coupling constants);
2. Hessians, or gradients, at geometries displaced along the normal modes
   (:func:`displacements`), whatever code or model computes them;
3. the cubic and semi-diagonal quartic force constants by finite differences
   (:func:`force_constants_from_hessians`, :func:`force_constants_from_gradients`);
4. VPT2 with Fermi resonances deperturbed, optionally treated variationally
   (GVPT2) (:func:`vpt2`).

Units: atomic units in, except masses in amu; the force constants and all
results in cm^-1, in dimensionless (reduced) normal coordinates
q_i = sqrt(omega_i) Q_i.

Only molecules without degenerate modes (asymmetric tops) are handled: the
formulas for degenerate modes, needed by linear molecules and symmetric and
spherical tops, are not yet in, so :class:`Harmonic` refuses those molecules.

The VPT2 formulas, the detection of Fermi resonances and the finite-difference
assembly are adapted from pyVPT2 (Copyright 2021-2024 Philip Nelson,
BSD-3-Clause; see THIRD_PARTY_LICENSES.rst).
"""

import itertools

import numpy as np

from .fermi import Interaction, State, fermi_solver

# CODATA 2018
HARTREE_TO_CM = 219474.6313632
HARTREE_TO_KJMOL = 2625.4996394799
AMU_TO_ME = 1822.888486209
BOHR_TO_ANGSTROM = 0.529177210903
CM_TO_KJMOL = HARTREE_TO_KJMOL / HARTREE_TO_CM

#: The steps, in units of the displacement, for each finite-difference route
STEPS = {
    "hessians": {"single": (-1, 1), "pair": ()},
    "gradients": {
        "single": (-3, -1, 1, 3),
        "pair": ((-1, -1), (1, -1), (-1, 1), (1, 1)),
    },
}


class DegenerateModesError(ValueError):
    """The molecule has degenerate vibrations, which VPT2 here cannot treat."""


class NotAMinimumError(ValueError):
    """The structure has imaginary or vanishing vibrational frequencies."""


class Harmonic:
    """The harmonic analysis at a minimum.

    Parameters
    ----------
    hessian : array (3N, 3N)
        The Cartesian Hessian, E_h/bohr^2.
    masses : array (N,)
        The atomic masses, amu.
    coordinates : array (N, 3)
        The coordinates, bohr.
    degenerate_tolerance : float
        Relative tolerance for equal principal moments of inertia, which mark a
        linear molecule or a symmetric or spherical top.

    Attributes
    ----------
    omega : array (n,)
        Harmonic frequencies, cm^-1, ascending.
    L : array (3N, n)
        Orthonormal mass-weighted normal modes.
    modes : array (3N, n)
        The Cartesian displacement, bohr, for a unit step in each dimensionless
        normal coordinate.
    B : array (3,)
        Rotational constants, cm^-1, A >= B >= C.
    zeta : array (3, n, n)
        Coriolis coupling constants about the principal axes of B.
    """

    def __init__(self, hessian, masses, coordinates, degenerate_tolerance=1.0e-3):
        H = np.asarray(hessian, dtype=float)
        masses = np.asarray(masses, dtype=float)
        xyz = np.asarray(coordinates, dtype=float).reshape(-1, 3)
        n_atoms = len(masses)
        if H.shape != (3 * n_atoms, 3 * n_atoms):
            raise ValueError(f"The Hessian is {H.shape}, not 3N x 3N for N={n_atoms}")
        H = 0.5 * (H + H.T)

        m = masses * AMU_TO_ME
        sqrt_m = np.repeat(np.sqrt(m), 3)
        center = (m[:, None] * xyz).sum(axis=0) / m.sum()
        r = xyz - center
        self.masses = masses
        self.coordinates = xyz

        # Principal moments and axes, m_e bohr^2
        inertia = np.zeros((3, 3))
        for mass, (x, y, z) in zip(m, r):
            inertia += mass * np.array(
                [
                    [y * y + z * z, -x * y, -x * z],
                    [-x * y, x * x + z * z, -y * z],
                    [-x * z, -y * z, x * x + y * y],
                ]
            )
        moments, axes = np.linalg.eigh(inertia)
        self.moments = moments
        self.axes = axes
        self.rotor = _rotor_type(moments, degenerate_tolerance)
        if self.rotor != "asymmetric top":
            raise DegenerateModesError(
                f"The molecule is a {self.rotor}, so it has degenerate vibrations. "
                "VPT2 here handles only asymmetric tops so far: the formulas for "
                "degenerate modes are not yet implemented."
            )
        self.B = HARTREE_TO_CM / (2 * moments)

        # Translations and rotations in mass-weighted coordinates, and the
        # orthonormal basis of the vibrations orthogonal to them
        tr = []
        for a in range(3):
            t = np.zeros((n_atoms, 3))
            t[:, a] = 1.0
            tr.append(t.ravel() * sqrt_m)
        for a in range(3):
            axis = axes[:, a]
            tr.append(np.cross(axis, r).ravel() * sqrt_m)
        T, _ = np.linalg.qr(np.array(tr).T)
        full, _, _ = np.linalg.svd(T, full_matrices=True)
        D = full[:, 6:]

        Hm = H / np.outer(sqrt_m, sqrt_m)
        values, U = np.linalg.eigh(D.T @ Hm @ D)
        if np.any(values <= 0) or np.any(np.sqrt(values) * HARTREE_TO_CM < 1.0):
            low = np.sign(values) * np.sqrt(np.abs(values)) * HARTREE_TO_CM
            raise NotAMinimumError(
                "The structure is not a minimum: lowest vibrational frequencies "
                + ", ".join(f"{w:.1f}" for w in low[:3])
                + " cm^-1 (negative = imaginary)."
            )
        # A deterministic phase for each mode (its largest component positive),
        # so the displacements are the same on every machine and on a restart
        L = D @ U
        for i in range(L.shape[1]):
            if L[np.argmax(np.abs(L[:, i])), i] < 0:
                L[:, i] = -L[:, i]
        self.L = L
        omega_au = np.sqrt(values)
        self.omega = omega_au * HARTREE_TO_CM
        self.modes = self.L / sqrt_m[:, None] / np.sqrt(omega_au)[None, :]
        self.n = len(self.omega)

        # Coriolis constants about the principal axes
        epsilon = np.zeros((3, 3, 3))
        epsilon[0, 1, 2] = epsilon[1, 2, 0] = epsilon[2, 0, 1] = 1.0
        epsilon[0, 2, 1] = epsilon[1, 0, 2] = epsilon[2, 1, 0] = -1.0
        self.zeta = np.zeros((3, self.n, self.n))
        for a in range(3):
            M = axes @ epsilon[a] @ axes.T
            self.zeta[a] = self.L.T @ np.kron(np.eye(n_atoms), M) @ self.L

    @property
    def zpe(self):
        """The harmonic zero-point energy, cm^-1."""
        return 0.5 * self.omega.sum()

    def to_normal(self, hessian):
        """A Cartesian Hessian (E_h/bohr^2) in dimensionless normal coordinates,
        cm^-1."""
        H = np.asarray(hessian, dtype=float)
        return self.modes.T @ (0.5 * (H + H.T)) @ self.modes * HARTREE_TO_CM

    def gradient_to_normal(self, gradient):
        """A Cartesian gradient (E_h/bohr) in dimensionless normal coordinates,
        cm^-1."""
        return self.modes.T @ np.asarray(gradient, dtype=float).ravel() * HARTREE_TO_CM


def _rotor_type(moments, tolerance):
    """'linear', 'spherical top', 'symmetric top' or 'asymmetric top' from the
    principal moments, ascending."""
    largest = moments[-1]
    if moments[0] < tolerance * largest:
        return "linear"
    equal = [
        abs(moments[i] - moments[j]) <= tolerance * largest for i, j in ((0, 1), (1, 2))
    ]
    if all(equal):
        return "spherical top"
    if any(equal):
        return "symmetric top"
    return "asymmetric top"


def displacements(harmonic, route="hessians", step=0.05):
    """The displaced geometries for the finite differences.

    Parameters
    ----------
    harmonic : Harmonic
    route : str
        'hessians' (2n geometries) or 'gradients' (2n^2 + 2n).
    step : float
        The displacement in dimensionless normal coordinates.

    Returns
    -------
    [(key, coordinates)]
        ``key`` is a tuple of (mode, steps) pairs, e.g. ``((3, 1),)`` or
        ``((1, -1), (4, 1))``; coordinates (N, 3) in bohr.
    """
    steps = STEPS[route]
    result = []
    for i in range(harmonic.n):
        for s in steps["single"]:
            result.append((((i, s),), _displace(harmonic, ((i, s),), step)))
    for i, j in itertools.combinations(range(harmonic.n), 2):
        for s, t in steps["pair"]:
            key = ((i, s), (j, t))
            result.append((key, _displace(harmonic, key, step)))
    return result


def _displace(harmonic, key, step):
    xyz = harmonic.coordinates.ravel().copy()
    for mode, steps in key:
        xyz += harmonic.modes[:, mode] * step * steps
    return xyz.reshape(-1, 3)


def force_constants_from_hessians(harmonic, hessians, step=0.05, reference=None):
    """The cubic and semi-diagonal quartic force constants from the Hessians at
    the displaced geometries.

    Parameters
    ----------
    harmonic : Harmonic
    hessians : dict
        Cartesian Hessians (E_h/bohr^2) by the keys from :func:`displacements`.
    step : float
        The displacement used.
    reference : array, optional
        The Hessian at the reference geometry at the same level as the displaced
        ones (mixed levels); by default the harmonic one's diagonal, omega.

    Returns
    -------
    phi_ijk : array (n, n, n)
    phi_iijj : array (n, n)
        In cm^-1, symmetrized.
    spread : dict
        The largest difference between finite-difference routes that should
        agree ('cubic', 'quartic'), cm^-1: a measure of the numerical noise.
    """
    n = harmonic.n
    h0 = harmonic.omega if reference is None else np.diag(harmonic.to_normal(reference))
    plus = {i: harmonic.to_normal(hessians[((i, 1),)]) for i in range(n)}
    minus = {i: harmonic.to_normal(hessians[((i, -1),)]) for i in range(n)}

    phi_ijk = np.zeros((n, n, n))
    phi_iijj = np.zeros((n, n))
    for i in range(n):
        # Displace along i, read the (j, k) element; the other permutations of
        # (i, j, k) are independent routes through j and k.
        phi_ijk[i] = (plus[i] - minus[i]) / (2 * step)
        # Displace along i, read the (j, j) element: phi_jjii
        phi_iijj[:, i] = (np.diag(plus[i]) + np.diag(minus[i]) - 2 * h0) / step**2
    return _symmetrize(phi_ijk, phi_iijj)


def force_constants_from_gradients(harmonic, gradients, step=0.05, reference=None):
    """The cubic and semi-diagonal quartic force constants from the gradients
    at the displaced geometries.

    Parameters
    ----------
    harmonic : Harmonic
    gradients : dict
        Cartesian gradients (E_h/bohr) by the keys from :func:`displacements`
        with route 'gradients'.
    step : float
        The displacement used.
    reference : array, optional
        The gradient at the reference geometry (zero at the minimum by default).

    Returns
    -------
    As :func:`force_constants_from_hessians`.
    """
    n = harmonic.n
    g = {key: harmonic.gradient_to_normal(value) for key, value in gradients.items()}
    g0 = np.zeros(n) if reference is None else harmonic.gradient_to_normal(reference)

    def get(*key):
        key = tuple(sorted(key))
        return g[key]

    phi_ijk = np.zeros((n, n, n))
    phi_iijj = np.zeros((n, n))
    for i in range(n):
        gp, gn = get((i, 1)), get((i, -1))
        g3p, g3n = get((i, 3)), get((i, -3))
        phi_iijj[i, i] = (g3p[i] - 3 * gp[i] + 3 * gn[i] - g3n[i]) / (8 * step**3)
        # Second difference of the gradient along i: phi_kii for every k, in
        # all three of its positions (there is only this one route to it)
        phi_kii = (gp + gn - 2 * g0) / step**2
        phi_ijk[:, i, i] = phi_kii
        phi_ijk[i, :, i] = phi_kii
        phi_ijk[i, i, :] = phi_kii
    for i, j in itertools.permutations(range(n), 2):
        pp, pn = get((i, 1), (j, 1)), get((i, 1), (j, -1))
        np_, nn = get((i, -1), (j, 1)), get((i, -1), (j, -1))
        gp, gn = get((i, 1)), get((i, -1))
        phi_iijj[i, j] = (pp[i] + pn[i] - 2 * gp[i] - np_[i] - nn[i] + 2 * gn[i]) / (
            2 * step**3
        )
    for i, j, k in itertools.permutations(range(n), 3):
        pp, nn = get((j, 1), (k, 1)), get((j, -1), (k, -1))
        pj, nj = get((j, 1)), get((j, -1))
        pk, nk = get((k, 1)), get((k, -1))
        phi_ijk[i, j, k] = (
            pp[i] + nn[i] + 2 * g0[i] - pj[i] - pk[i] - nj[i] - nk[i]
        ) / (2 * step**2)
    return _symmetrize(phi_ijk, phi_iijj)


def _symmetrize(phi_ijk, phi_iijj):
    """Average the routes that should agree, and report how far apart they
    were."""
    cubic = max(
        (
            np.max(np.abs(phi_ijk - phi_ijk.transpose(p)))
            for p in itertools.permutations(range(3))
        ),
        default=0.0,
    )
    quartic = float(np.max(np.abs(phi_iijj - phi_iijj.T))) if phi_iijj.size else 0.0
    sym = sum(phi_ijk.transpose(p) for p in itertools.permutations(range(3))) / 6
    return (
        sym,
        0.5 * (phi_iijj + phi_iijj.T),
        {"cubic": float(cubic), "quartic": quartic},
    )


def find_resonances(omega, phi_ijk, omega_threshold=200.0, k_threshold=1.0):
    """The Fermi resonances, by Martin's criteria.

    Returns
    -------
    [dict]
        Each with 'type' (1: omega_j ~ 2 omega_i; 2: omega_k ~ omega_i +
        omega_j), 'fundamental' (the mode), 'states' (the pair of modes),
        'delta' (the mismatch, cm^-1) and 'K' (the second-order shift, cm^-1).
    """
    n = len(omega)
    found = []
    for i, j in itertools.permutations(range(n), 2):
        d = abs(2 * omega[i] - omega[j])
        if d <= omega_threshold:
            K = phi_ijk[i, i, j] ** 4 / (256 * d**3) if d > 0 else np.inf
            if K >= k_threshold:
                found.append(
                    {"type": 1, "fundamental": j, "states": (i, i), "delta": d, "K": K}
                )
    for i, j, k in itertools.permutations(range(n), 3):
        if i > j:
            continue
        d = abs(omega[i] + omega[j] - omega[k])
        if d <= omega_threshold:
            K = phi_ijk[i, j, k] ** 4 / (64 * d**3) if d > 0 else np.inf
            if K >= k_threshold:
                found.append(
                    {"type": 2, "fundamental": k, "states": (i, j), "delta": d, "K": K}
                )
    return found


def vpt2(
    harmonic,
    phi_ijk,
    phi_iijj,
    fermi=True,
    gvpt2=True,
    omega_threshold=200.0,
    k_threshold=1.0,
    low_mode=200.0,
):
    """VPT2 frequencies and zero-point energy.

    Parameters
    ----------
    harmonic : Harmonic
    phi_ijk, phi_iijj : arrays
        The cubic and semi-diagonal quartic force constants, cm^-1.
    fermi : bool
        Deperturb Fermi resonances (Martin's criteria).
    gvpt2 : bool
        Then treat them variationally (GVPT2). Needs ``fermi``.
    omega_threshold, k_threshold : float
        Martin's criteria for a resonance, cm^-1.
    low_mode : float
        Modes below this, cm^-1, are flagged: VPT2 is unreliable for torsions
        and other large-amplitude motions.

    Returns
    -------
    dict
        'omega', 'nu' (fundamentals), 'nu_deperturbed', 'overtones',
        'combinations' (n x n), 'chi' (n x n), 'zpe_harmonic', 'zpe' (the
        anharmonic ZPE including the constant term), 'g0' (the constant term),
        'resonances', 'flagged' (mode indices); all cm^-1.
    """
    omega = harmonic.omega
    n = len(omega)
    B, zeta = harmonic.B, harmonic.zeta
    resonances = (
        find_resonances(omega, phi_ijk, omega_threshold, k_threshold) if fermi else []
    )
    type1 = {(r["fundamental"], r["states"][0]) for r in resonances if r["type"] == 1}
    type2 = set()
    for r in resonances:
        if r["type"] == 2:
            k, (i, j) = r["fundamental"], r["states"]
            type2.add((k, i, j))
            type2.add((k, j, i))

    def resonant(f, a, b):
        """Whether omega_f ~ omega_a + omega_b is a treated resonance: type 2,
        or type 1 when a == b (omega_f ~ 2 omega_a), which enters the
        off-diagonal chi through its k = i and k = j terms."""
        return (f, a, b) in type2 or (a == b and (f, a) in type1)

    chi = np.zeros((n, n))
    chi0 = 0.0
    for i in range(n):
        chi0 += phi_iijj[i, i]
        chi0 -= (7 / 9) * phi_ijk[i, i, i] ** 2 / omega[i]

        # Diagonal
        value = phi_iijj[i, i]
        for k in range(n):
            if (k, i) in type1:  # omega_k ~ 2 omega_i
                value -= (phi_ijk[i, i, k] ** 2 / 2) * (
                    1 / (2 * omega[i] + omega[k]) + 4 / omega[k]
                )
            else:
                value -= (
                    (8 * omega[i] ** 2 - 3 * omega[k] ** 2)
                    * phi_ijk[i, i, k] ** 2
                    / (omega[k] * (4 * omega[i] ** 2 - omega[k] ** 2))
                )
        chi[i, i] = value / 16

        # Off-diagonal
        for j in range(n):
            if j == i:
                continue
            if (i, j) in type1:  # omega_i ~ 2 omega_j
                chi0 -= (
                    3
                    * omega[j]
                    * phi_ijk[i, j, j] ** 2
                    / (4 * omega[i] ** 2 + 2 * omega[i] * omega[j])
                )
            else:
                chi0 += (
                    3
                    * omega[i]
                    * phi_ijk[i, j, j] ** 2
                    / (4 * omega[j] ** 2 - omega[i] ** 2)
                )
            value = phi_iijj[i, j]
            rot = sum(B[b] * zeta[b, i, j] ** 2 for b in range(3))
            value += 4 * (omega[i] ** 2 + omega[j] ** 2) / (omega[i] * omega[j]) * rot
            for k in range(n):
                value -= phi_ijk[i, i, k] * phi_ijk[j, j, k] / omega[k]
                wi, wj, wk = omega[i], omega[j], omega[k]
                if resonant(k, i, j):  # omega_k ~ omega_i + omega_j
                    d_ij = (
                        1 / (wi + wj + wk) + 1 / (-wi + wj + wk) + 1 / (wi - wj + wk)
                    ) / -2
                    d_0 = 1 / (wi + wj + wk) - 1 / (wi - wj + wk) - 1 / (-wi + wj + wk)
                elif resonant(i, j, k):  # omega_i ~ omega_j + omega_k
                    d_ij = (
                        1 / (wi + wj + wk) - 1 / (wi + wj - wk) + 1 / (wi - wj + wk)
                    ) / -2
                    d_0 = 1 / (wi + wj + wk) - 1 / (wi + wj - wk) - 1 / (wi - wj + wk)
                elif resonant(j, i, k):  # omega_j ~ omega_i + omega_k
                    d_ij = (
                        1 / (wi + wj + wk) - 1 / (wi + wj - wk) + 1 / (-wi + wj + wk)
                    ) / -2
                    d_0 = 1 / (wi + wj + wk) - 1 / (wi + wj - wk) - 1 / (-wi + wj + wk)
                else:
                    delta = (
                        (wi + wj - wk)
                        * (wi + wj + wk)
                        * (wi - wj + wk)
                        * (wi - wj - wk)
                    )
                    d_ij = 2 * wk * (wi**2 + wj**2 - wk**2) / delta
                    d_0 = -8 * wi * wj * wk / delta
                value += phi_ijk[i, j, k] ** 2 * d_ij
                if j > i and k > j:
                    chi0 += 2 * phi_ijk[i, j, k] ** 2 * d_0
            chi[i, j] = value / 4

    for b in range(3):
        zeta_sum = sum(
            zeta[b, i, j] ** 2 for i, j in itertools.combinations(range(n), 2)
        )
        chi0 -= 16 * B[b] * (1 + 2 * zeta_sum)
    g0 = chi0 / 64

    zpe = g0
    for i in range(n):
        zpe += 0.5 * (omega[i] + 0.5 * chi[i, i])
        for j in range(i + 1, n):
            zpe += 0.25 * chi[i, j]

    nu = np.zeros(n)
    overtones = np.zeros(n)
    for i in range(n):
        nu[i] = omega[i] + 2 * chi[i, i]
        overtones[i] = 2 * omega[i] + 6 * chi[i, i]
        for j in range(n):
            if j != i:
                nu[i] += 0.5 * chi[i, j]
                overtones[i] += chi[i, j]
    combinations = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        value = omega[i] + omega[j] + 2 * chi[i, i] + 2 * chi[j, j] + 2 * chi[i, j]
        for k in range(n):
            if k not in (i, j):
                value += 0.5 * (chi[i, k] + chi[j, k])
        combinations[i, j] = combinations[j, i] = value

    deperturbed = nu.copy()
    if fermi and gvpt2 and resonances:
        interactions = []
        for r in resonances:
            k = r["fundamental"]
            i, j = r["states"]
            if r["type"] == 1:
                right = State(state=(i, i), nu=overtones[i])
            else:
                right = State(state=(i, j), nu=combinations[i, j])
            interactions.append(
                Interaction(
                    left=State(state=(k,), nu=nu[k]),
                    right=right,
                    phi=phi_ijk[k, i, j],
                    ftype=r["type"],
                )
            )
        for state, value in fermi_solver(interactions).items():
            if len(state) == 1:
                nu[state[0]] = value

    return {
        "omega": omega,
        "nu": nu,
        "nu_deperturbed": deperturbed,
        "overtones": overtones,
        "combinations": combinations,
        "chi": chi,
        "zpe_harmonic": harmonic.zpe,
        "zpe": zpe,
        "g0": g0,
        "resonances": resonances,
        "flagged": [i for i in range(n) if omega[i] < low_mode],
    }

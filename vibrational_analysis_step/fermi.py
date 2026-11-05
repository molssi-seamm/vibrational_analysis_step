# -*- coding: utf-8 -*-

"""The variational treatment of Fermi resonances in generalized VPT2 (GVPT2).

The deperturbed VPT2 energies of the states in a resonance -- a fundamental and
the overtone or combination band it is resonant with -- are the diagonal of a
small effective Hamiltonian, and the cubic force constant that couples them its
off-diagonal elements. States joined by any chain of resonances form one polyad,
whose Hamiltonian is diagonalized; each state takes the eigenvalue of the
eigenvector it dominates.

Adapted from pyVPT2's ``fermi_solver.py`` (Copyright 2021-2024 Philip Nelson,
BSD-3-Clause; see THIRD_PARTY_LICENSES.rst).
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class State:
    """A vibrational state and its (deperturbed) frequency.

    Attributes
    ----------
    state : tuple of int
        The modes excited: ``(i,)`` a fundamental, ``(j, j)`` an overtone,
        ``(j, k)`` a combination band.
    nu : float
        Its frequency in cm^-1.
    """

    state: tuple
    nu: float


@dataclass
class Interaction:
    """A Fermi resonance between two states.

    Attributes
    ----------
    left, right : State
        The fundamental, and the overtone or combination band.
    phi : float
        The cubic force constant coupling them, in cm^-1.
    ftype : int
        1 for a fundamental and an overtone (i = 2j), 2 for a fundamental and a
        combination band (i = j + k).
    """

    left: State
    right: State
    phi: float
    ftype: int


def fermi_solver(interactions):
    """The variationally corrected frequencies of the states in resonance.

    Parameters
    ----------
    interactions : [Interaction]

    Returns
    -------
    dict
        Each state in a resonance (a tuple of mode indices) and its frequency.
    """
    polyads = []
    for interaction in interactions:
        # Every polyad already holding one of the two states; an interaction
        # that bridges several merges them.
        hits = [
            p
            for p in polyads
            if interaction.left.state in p.states or interaction.right.state in p.states
        ]
        if not hits:
            polyads.append(Polyad(interaction))
            continue
        target = hits[0]
        for other in hits[1:]:
            target.merge(other)
            polyads.remove(other)
        target.add(interaction)

    result = {}
    for polyad in polyads:
        result.update(polyad.solve())
    return result


class Polyad:
    """States coupled by Fermi resonances, and their effective Hamiltonian."""

    def __init__(self, interaction):
        self.states = set()
        self.nu = {}
        self.couplings = {}
        self.add(interaction)

    def add(self, interaction):
        """Add an interaction and its two states."""
        left, right = interaction.left, interaction.right
        self.states.update([left.state, right.state])
        self.nu[left.state] = left.nu
        self.nu[right.state] = right.nu
        self.couplings[(left.state, right.state)] = (interaction.phi, interaction.ftype)

    def merge(self, other):
        """Take over the states and interactions of another polyad."""
        self.states.update(other.states)
        self.nu.update(other.nu)
        self.couplings.update(other.couplings)

    def hamiltonian(self):
        """The effective Hamiltonian and the index of each state in it."""
        index = {state: i for i, state in enumerate(sorted(self.states))}
        H = np.zeros((len(index), len(index)))
        for state, i in index.items():
            H[i, i] = self.nu[state]
        for (a, b), (phi, ftype) in self.couplings.items():
            i, j = index[a], index[b]
            if ftype == 1:
                H[i, j] = H[j, i] = phi / 4
            else:
                H[i, j] = H[j, i] = phi / (2 * np.sqrt(2))
        return H, index

    def solve(self):
        """Diagonalize; each state gets the eigenvalue of the eigenvector it
        dominates.

        One eigenvalue per state: repeatedly take the largest remaining squared
        coefficient (state, eigenvector) and retire both, so no two states can
        claim the same eigenvector.
        """
        H, index = self.hamiltonian()
        values, vectors = np.linalg.eigh(H)
        weight = np.square(vectors)
        of_state = {}
        for _ in range(len(values)):
            row, column = np.unravel_index(np.argmax(weight), weight.shape)
            of_state[row] = values[column]
            weight[row, :] = -1.0
            weight[:, column] = -1.0
        return {state: of_state[i] for state, i in index.items()}

.. _user-guide:

**********
User Guide
**********

The Vibrational Analysis step computes **anharmonic vibrational frequencies and the
anharmonic zero-point energy** of a molecule by second-order vibrational perturbation
theory (VPT2), with whatever method a **Model Chemistry** step earlier in the flowchart
chooses: ORCA, Gaussian, MOPAC, xTB, an MLFF, ... It does the VPT2 itself rather than
relying on each program's own implementation, so it works, and gives the same kind of
results, with every program SEAMM drives.

What it gives
=============

Harmonic frequencies overestimate both the zero-point energy and the observed
fundamentals, and no single scale factor corrects both, or corrects different kinds of
molecules consistently. VPT2 adds the cubic and quartic force constants:

* the **anharmonic fundamentals**, the **overtones** (2ν\ :sub:`i`) and
  **combination bands** (ν\ :sub:`i` + ν\ :sub:`j`), comparable with experiment band for
  band;
* the **anharmonic zero-point energy**, including the constant term G\ :sub:`0`, the
  same convention as Gaussian's;
* the anharmonicity constants χ.

**Limitations.** Molecules with degenerate vibrations -- linear molecules and symmetric
or spherical tops such as CH\ :sub:`4`, NH\ :sub:`3`, ethane or benzene -- are refused
with a clear message: the formulas for degenerate modes are not in yet. The structure
must be a minimum. Modes below a threshold (200 cm\ :sup:`-1` by default), torsions and
other large-amplitude motions, are flagged in the output: VPT2 is not reliable for them,
though they are kept in the zero-point energy.

How to use it
=============

Put a **Model Chemistry** step before it, and the molecule in the current structure:

.. code-block:: yaml

   steps:
   - Model Chemistry: {model chemistry: "ORCA:DFT@B3LYP/def2-TZVP"}
   - FromSMILESStep: {smiles string: O}
   - Vibrational Analysis

By default the step optimizes the structure tightly first, over the program's warm MDI
engine, then computes:

1. the Hessian at the minimum, giving the normal modes and harmonic frequencies;
2. Hessians, or gradients, at structures displaced along the normal modes;
3. the cubic and semi-diagonal quartic force constants, by finite differences;
4. VPT2, with Fermi resonances found by Martin's criteria, removed from the
   perturbation sums, and the resonant states mixed variationally (GVPT2).

Hessians or gradients
---------------------

**At the displaced structures** chooses what is computed there. ``automatic`` takes
Hessians when the program has analytic second derivatives (ORCA for HF, most DFT
functionals and MP2), which needs 2n Hessians for n modes (54 for propane), and
gradients otherwise, about 2n\ :sup:`2` of them (about 1,500 for propane), each much
cheaper. Both give the same results.

**Calculations via** chooses how the program runs: ``automatic`` uses the warm MDI
engine when the program has one and is installed here, else separate calculations
where the job runs them; ``separate tasks`` runs them in parallel (on this machine or a
cluster queue) but cannot optimize the structure, so optimize it with the program's
own step first.

Precision and step size
-----------------------

Finite differences magnify errors in the Hessians or gradients: the quartic constants
are second differences of Hessians, or third differences of gradients. The step
therefore asks the program for **tight convergence** (xTB accuracy 1e-4, MOPAC relative
tolerance 0.001, ORCA VeryTightSCF), both for the MDI engine and for separate
calculations, and reports the **numerical noise**: how far apart the finite-difference
routes to the same force constant are. If the quartic noise exceeds 5 cm\ :sup:`-1` it
warns that the results may be off by several cm\ :sup:`-1`.

**Displacement** is the step along each mode, in dimensionless normal coordinates.
``automatic`` is 0.20 for Hessians and 0.10 for gradients; ``mass-weighted`` is
Gaussian's convention, 0.025 Å amu\ :sup:`1/2` for every mode. Smaller steps magnify
the noise, larger ones add higher-order error.

Resonances and low modes
------------------------

**Fermi resonances**: ``GVPT2 (variational)`` (the default), ``deperturbed VPT2``, or
``none (plain VPT2)``, which can give wildly wrong frequencies near a resonance. The
**resonance window** (200 cm\ :sup:`-1`) and **threshold** (1 cm\ :sup:`-1`) are
Martin's two criteria. The zero-point energy does not depend on this choice.

Results
=======

The output lists, for each mode, the harmonic and anharmonic frequency, the shift and the
overtone, marking resonances and low-frequency modes, then the zero-point energies and
the noise. The force constants and all results are written to ``vpt2.json`` in the
step's directory. Results that can be saved to variables, tables or properties:

* ``ZPE`` (property ``ZPE#VPT2#{model}``), ``harmonic ZPE``, ``ZPE correction`` (kJ/mol);
* ``fundamentals``, ``harmonic frequencies``, ``overtones``, ``combination bands`` and
  the anharmonicity constants ``chi`` (cm\ :sup:`-1`, json properties);
* ``G0``, the resonances treated, the flagged modes and the number of modes.

Validation
==========

Against Gaussian 09's own anharmonic calculation (B3LYP/6-31G(d), Gaussian Hessians at
this step's displacements), the zero-point energies agree to 0.3 cm\ :sup:`-1` for water
and propane and 2 cm\ :sup:`-1` for formaldehyde, the fundamentals to 0.7 cm\ :sup:`-1`
rms for propane. Against ORCA 6.1.1's own VPT2 for water (B3LYP/def2-TZVP) the
fundamentals agree to 0.02 cm\ :sup:`-1`.

References
==========

The VPT2 formulas are those of Mills (1972); the resonance criteria are Martin et al.'s
(1995) and the generalized treatment Barone's (2005), with the polyad handling of Franke
et al. (2021). The algebra is adapted from pyVPT2 (Philip Nelson, BSD-3-Clause). All are
cited in the job's references.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

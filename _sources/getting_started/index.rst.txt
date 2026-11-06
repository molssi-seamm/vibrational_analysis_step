***************
Getting Started
***************

Installation
============
The Vibrational Analysis step is installed with SEAMM's other plug-ins. To add or
update it in an installation::

  seamm-manager install vibrational-analysis-step

An example
==========
The anharmonic frequencies and zero-point energy of water with ORCA, in a flowchart
built from a spec (see the SEAMM documentation on building flowcharts):

.. code-block:: yaml

   title: Water VPT2
   steps:
   - Model Chemistry: {model chemistry: "ORCA:DFT@B3LYP/def2-SVP"}
   - FromSMILESStep: {smiles string: O}
   - Vibrational Analysis

The step optimizes the structure, computes the analytic Hessian and the Hessians at six
displaced structures over ORCA's MDI engine, and prints the fundamentals (here about
1583, 3609 and 3698 cm\ :sup:`-1`, against harmonic 1638, 3788 and 3883) and the
anharmonic zero-point energy (a correction of about -74 cm\ :sup:`-1`). See the
:ref:`User Guide <user-guide>` for the options.

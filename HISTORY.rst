=======
History
=======

2026.11.5.1 -- Zero-point energies of symmetric and spherical tops
    * Molecules with degenerate vibrations, symmetric tops such as NH3 and ethane and
      spherical tops such as CH4, are no longer refused: the step gives their
      anharmonic zero-point energy and harmonic frequencies. The anharmonic fundamentals
      and constants of degenerate modes are not available yet. Linear molecules are
      still refused.
    * Validated against Gaussian 09 to about 1 cm^-1 by giving the molecules
      slightly asymmetric masses, where Gaussian uses the general formulas.

2026.11.5 -- Initial release: anharmonic frequencies (VPT2) with any Model Chemistry
    * Anharmonic fundamentals, overtones, combination bands and the anharmonic zero-point
      energy by second-order vibrational perturbation theory, with Fermi resonances
      treated variationally (GVPT2), for any program a Model Chemistry step chooses.
    * Hessians or gradients at structures displaced along the normal modes, over the
      program's warm MDI engine or as separate calculations in parallel, with tight
      convergence and a warning when the force constants are noisy.
    * Validated against Gaussian 09's and ORCA 6.1.1's own VPT2. Molecules with
      degenerate vibrations (linear molecules, symmetric and spherical tops) are not
      handled yet.

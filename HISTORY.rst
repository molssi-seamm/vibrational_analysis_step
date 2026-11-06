=======
History
=======

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

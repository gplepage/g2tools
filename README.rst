g2tools
-------

This module contains a small number of tools useful for analyzing
contributions to the muon's magnetic moment from (lattice) QCD vacuum
polarization. The functions or classes include:

    **moments(G)**
        Compute moments of *jj* correlator *G*.

    **mom2taylor(mom)**
        Convert moments ``mom`` into Taylor coefficients for *q*\ :sup:`2`-expansion.

    **taylor2mom(tayl)**
        Convert Taylor expansion coefficients ``tayl`` into moments.

    **vacpol(mom)**
        Create a Pade approximant for the subtracted
        vacuum polarization (``PI-hat``) from the *jj* correlator
        whose moments (or Taylor coefficients) are in *mom*.

    **fourier_vacpol(G)**
        Create subtracted vacuum polarization (``PI-hat``) by
        Fourier transforming *jj* correlator ``G(t)``.

    **a_mu(pihat, Q)**
        Compute the contribution to the muon's *g-2*
        anomaly from function pihat (usually built by vacpol).

    **R2G(E, R)**
        Compute the Euclidean G(t) corresponding to data 
        for Re+e-.

    **R2a_mu(E, R)**
        Compute the leading-order contribution to the 
        muon's *g-2* anomaly corresponding to data 
        for Re+e-.

    **TanhWin(t0, t1, dt)**
        Create a filter for applying a t-window in 
        ``monents(...)`` or ``fourier_vacpol(...)``.

    **pade_gvar(f, m, n)**
        General-purpose code for determining Pade approximants
        to a power series whose coefficients are ``GVar``\s (ie,
        Gaussian random variables, for error propagation).

    **pade_svd(f, m, n)**
        General-purpose code for determining Pade approximants
        for a power series whose coefficients are floats.
        Uses *svd* regularization to stabilize results when
        the input data are noisy.

Information on how to install the module is in the file INSTALLATION.

To test the module try ``make tests``.

Documentation is in the doc directory: open doc/html/index.html
or look online at <https://g2tools.readthedocs.io>.

The examples directory has a complete example, showing how to go from Monte
Carlo data for a *jj* correlator to a contribution to the muon's magnetic
moment anomaly *a*\ :sub:`µ`. See also the introduction in the documentation.

The general technique that underlies this module is described in
Chakraborty *et al*, Phys.Rev. D89 (2014) no.11, 114501. Google
``arXiv:1403.1778`` to find a preprint on the web.

| Created by G. Peter Lepage (Cornell University) on on 2014-09-13.
| Copyright (c) 20014-22 G. Peter Lepage.

.. image:: https://zenodo.org/badge/66222496.svg
   :target: https://zenodo.org/badge/latestdoi/66222496


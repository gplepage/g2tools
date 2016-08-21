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

    **a_mu(pihat, Q)**
        Compute the contribution to the muon's *g-2*
        anomaly from function pihat (usually built by vacpol).

    **pade_gvar(f,m,n)**
        General-purpose code for determining Pade approximants
        to a power series whose coefficients are ``GVar``\s (ie,
        Gaussian random variables -- for error propagation).

    **pade_svd(f,m,n)**
        General-purpose code for determining Pade approximants
        for a power series whose coefficients are floats.
        Uses *svd* regularization to stabilize results when
        the input data are noisy.

Information on how to install the module is in the file INSTALLATION.

To test the module try ``make tests``.

Documentation is in the doc directory --- open doc/html/index.html.
A pdf version is in doc/g2tools.pdf.

The examples directory has a complete example, showing how to go from Monte
Carlo data for a *jj* correlator to a contribution to the muon's magnetic
moment anomaly *a*\ :sub:`µ`. See also the introduction in the documentation.

The general technique that underlies this module is described in
Chakraborty *et al*, Phys.Rev. D89 (2014) no.11, 114501. Google
``arXiv:1403.1778`` to find a preprint on the web.

Created by G. Peter Lepage (Cornell University) on on 2014-09-13.
Copyright (c) 20014-16 G. Peter Lepage.

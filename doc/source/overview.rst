Muon g-2 from Lattice QCD using :mod:`g2tools`
===============================================

.. |g2tools| replace:: :mod:`g2tools`

.. |GVar| replace:: :class:`gvar.GVar`

Module |g2tools| contains a small number of tools useful for analyzing
contributions to the muon's magnetic moment from (lattice) QCD vacuum
polarization. These tools were developed by G.P. Lepage to implement the
analysis presented in Chakraborty *et al*, Phys.Rev. D89 (2014) no.11, 114501
(arXiv:1403.1778) and subsequent papers by the same authors.

A typical application, illustrating the most important tools, is
provided by the following code::

    import g2tools as g2
    import gvar as gv

    def main():
        # data
        Z = gv.gvar('0.9938(17)')           # current Z factor
        Q = 1. / 3.                         # charge of quark (units of proton charge)
        ainv = gv.gvar('1.6280(86)')        # inverse lattice spacing (in GeV)

        G = gv.gvar([                       # G(t) for t=0..63 (in lattice units)
            '0.0870904(11)', '0.0435138(14)', '0.00509859(48)', '0.00305614(43)',
            '0.00069516(19)', '0.00045466(15)', '0.000166972(80)', '0.000102219(58)',
            '0.000045284(34)', '0.000026213(22)', '0.000012630(14)', '7.0635(91)e-06',
            '3.5569(57)e-06', '1.9469(37)e-06', '1.0027(24)e-06', '5.421(16)e-07',
            '2.834(10)e-07', '1.5174(67)e-07', '7.943(43)e-08', '4.253(28)e-08',
            '2.221(19)e-08', '1.183(12)e-08', '6.132(81)e-09', '3.292(51)e-09',
            '1.727(34)e-09', '9.19(22)e-10', '4.81(14)e-10', '2.643(96)e-10',
            '1.385(64)e-10', '7.61(44)e-11', '3.92(31)e-11', '2.67(24)e-11',
            '2.07(21)e-11', '2.90(23)e-11', '4.12(31)e-11', '8.20(42)e-11',
            '1.380(65)e-10', '2.788(98)e-10', '5.01(15)e-10', '9.72(23)e-10',
            '1.782(34)e-09', '3.406(53)e-09', '6.333(78)e-09', '1.212(12)e-08',
            '2.249(18)e-08', '4.283(28)e-08', '8.016(44)e-08', '1.5263(67)e-07',
            '2.843(10)e-07', '5.420(16)e-07', '1.0062(25)e-06', '1.9453(39)e-06',
            '3.5611(58)e-06', '7.0675(93)e-06', '0.000012647(14)', '0.000026240(22)',
            '0.000045282(32)', '0.000102285(56)', '0.000166993(79)', '0.00045479(15)',
            '0.00069503(19)', '0.00305647(42)', '0.00509870(47)', '0.0435158(14)'
            ])
        # N.B.: In general would construct G so that correlations from one t
        #   to the next are included. Don't bother here since this is meant
        #   just to illustrate g2tools.

        # compute moments, converting to physical units from lattice units
        mom = g2.moments(G, ainv=ainv, Z=Z, periodic=True, nlist=[4, 6, 8, 10])
        print('Taylor coefficients:', g2.mom2taylor(mom))
        print()

        # construct subtracted vac pol function using [2,2] Pade
        vpol = g2.vacpol(mom, order=(2,2))

        # integrate vpol to get a_mu and print result
        amu = g2.a_mu(vpol, Q=Q)
        print('a_mu contribution =', amu)
        print()

        # error budget for a_mu
        print(gv.fmt_errorbudget(
            outputs=dict(a_mu=a_mu, mom4=mom[4]),
            inputs=dict(G=G, Z=Z, ainv=ainv),
            ))

    if __name__ == '__main__':
        main()

In this code, we first read the simulation data for the *jj* correlator into
array ``G``, where ``G[i]`` is the correlator for (Euclidean) time  ``i/ainv``
where ``i=0,1..63``. We then use :func:`g2tools.moments` to calculate
temporal moments of the correlator, while also converting from lattice units
to physical units (using the inverse lattice spacing ``ainv``)  and
renormalizing the current (``Z``).

``vpol(q2)`` is the vacuum polarization function at Euclidean *q*\ :sup:`2`
equal to ``q2``. Object ``vpol`` has type :class:`g2tools.vacpol`. It
constructs a  [2,2] Padé approximant from the moments,
and uses that approximant to  approximate the exact function.
The approximants converge to the exact result as the order
increases provided the momentum is space-like (``q2`` non-negative).
Using a [1,1] Padé instead of [2,2] gives almost identical results here, so the
approximants have converged for the present application.

We calculate the contribution from vacuum polarization ``vpol``
to the muon's anomalous magnetic moment a\ :sub:`µ` using
:func:`g2tools.a_mu`. We also use :func:`gvar.fmt_errorbudget`
to produce an error budget for it and the 4th moment.

Running this code gives the following output::

    Taylor coefficients: [0.06629(74) -0.0527(11) 0.0472(15) -0.0435(18)]

    a_mu contribution = 5.412(57)e-09

    Partial % Errors:
                    a_mu      mom4
    ------------------------------
         ainv:      1.00      1.06
            Z:      0.34      0.34
            G:      0.01      0.01
    ------------------------------
        total:      1.06      1.11

The contribution to the muon's anomalous magnetic moment is
54.12(57)x10\ :sup:`-10`. The error budget shows that the final
uncertainty is dominated by the uncertainty in the inverse
lattice spacing ``ainv``; statistical errors from ``G`` are
completely negligible in this example.

|g2tools| is designed to work with module :mod:`gvar` which we use here
to represent the statistical and systematic uncertainties in
the correlator values, inverse lattice spacing, and ``Z`` factor. Each of these
quantities is an object of type |GVar|, which represents
a Gaussian random variable. |GVar|\s describe not only
means and standard deviations, but also statistical correlations between
different objects. These correlations are propagated through arbitrary
arithmetic statements. Adding the following code to the end of ``main()``,
for example, ::

    print(gv.evalcorr([mom[4], mom[6], mom[8], mom[10]]))

prints out the correlation matrix for the moments, showing that they
are highly correlated (as expected)::

    [[ 1.          0.98833867  0.9787737   0.97262094]
     [ 0.98833867  1.          0.99853653  0.99646438]
     [ 0.9787737   0.99853653  1.          0.99949934]
     [ 0.97262094  0.99646438  0.99949934  1.        ]]

The moments are also highly correlated with the final results ``a_mu``: for
example, adding the following to the end of ``main()`` ::

    print(gv.evalcorr([a_mu, mom[4]]))

gives::

    [[ 1.          0.96864247]
     [ 0.96864247  1.        ]]

This kind of correlation information is used by ``gvar.fmt_errorbudget(...)``
to create the error budget. See :mod:`gvar`'s documentation
for more information.


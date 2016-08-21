""" Example application of g2tools.py

This example illustrates the use of :mod:`g2tools` to convert
lattices data (adapted from Chakraborty et al, arXiv:1403.1778)
into a contribution to the muon's magnetic anomaly.
"""
from __future__ import print_function

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
    a_mu = g2.a_mu(vpol, Q=Q)
    print('a_mu contribution =', a_mu)
    print()

    # error budget for a_mu
    print(gv.fmt_errorbudget(
        outputs=dict(a_mu=a_mu, mom4=mom[4]),
        inputs=dict(G=G, Z=Z, ainv=ainv),
        ))

if __name__ == '__main__':
    main()


# Created by G. Peter Lepage (Cornell University) on 2016-08-20.
# Copyright (c) 2016 G. Peter Lepage.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version (see <http://www.gnu.org/licenses/>).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

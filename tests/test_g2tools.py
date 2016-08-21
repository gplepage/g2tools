"""
test_g2tools.py
"""
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

from __future__ import print_function   # makes this work for python2 and 3

import unittest
import numpy as np
import gvar as gv
from g2tools import *

SHOW_OUTPUT = False

def optprint(*args):
    pass

if SHOW_OUTPUT:
    optprint = print

MPI = 0.13957
MK = 0.4937

class test_g2tools(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_moments(self):
        " moments(G) "
        optprint('\n=========== Test moments')
        mom = moments([1., 2., 3., 2.], nlist=[0, 2, 4])
        assert mom[0] == 8. and mom[2] == 16. and mom[4] == 52.
        mom = moments([1., 2., 3., 2.], nlist=[0, 2, 4], periodic=False)
        assert mom[0] == 16. and mom[2] == 64. and mom[4] == 424.
        tayl = [1., -2. , 3.]
        mom = taylor2mom(tayl)
        assert mom[4] == 24. and mom[6] == 1440. and mom[8] == 120960.
        assert numpy.allclose(mom2taylor(mom), tayl)
        optprint('nothing to report -- all is good')

    def test_pade_svd(self):
        " pade_svd(tayl, n, m) "
        optprint('\n=========== Test pade_svd')
        from scipy.misc import pade
        # Taylor expansion for exp(x)
        e_exp = [1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0]

        # test against scipy
        p0, q0 = pade(e_exp, 2)
        p0 = p0.c[-1::-1]
        q0 = q0.c[-1::-1]
        p, q = pade_svd(e_exp, 3, 2)
        assert numpy.allclose(p, p0)
        assert numpy.allclose(q, q0)
        optprint('(3,2) Pade of exp(x) - num:', p)
        optprint('(3,2) Pade of exp(x) - den:', q)
        e = sum(p) / sum(q)
        optprint('Pade(x=1) = {:.6}    error = {:7.2}'.format(
            e,
            abs(e/numpy.exp(1) - 1.),
            ))

        # now with 10% errors --- automatically reduces to (2,1)
        p0, q0 = pade(e_exp[:4], 1)
        p0 = p0.c[-1::-1]
        q0 = q0.c[-1::-1]
        p, q = pade_svd(e_exp, 3, 2, rtol=0.1)
        assert numpy.allclose(p, p0)
        assert numpy.allclose(q, q0)
        optprint('(2,1) Pade of exp(x) - num:', p)
        optprint('(2,1) Pade of exp(x) - den:', q)
        e = sum(p) / sum(q)
        optprint('Pade(x=1) = {:.6}    error = {:7.2}'.format(
            e,
            abs(e/numpy.exp(1) - 1.)
            ))

        # now with 90% errors --- automatically reduces to (1,0)
        p, q = pade_svd(e_exp, 3, 2, rtol=0.9)
        optprint('(1,0) Pade of exp(x) - num:', p)
        optprint('(1,0) Pade of exp(x) - den:', q)
        e = sum(p) / sum(q)
        optprint('Pade(x=1) = {:.6}    error = {:7.2}'.format(
            e,
            abs(e/numpy.exp(1) - 1.)
            ))
        assert numpy.allclose(p, [1., 1.])
        assert numpy.allclose(q, [1.])

    def test_pade_gvar(self):
        " pade_gvar(tayl, m, n) "
        optprint('\n=========== Test pade_gvar')
        from scipy.misc import pade
        e_exp = [1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0]
        def scipy_pade(m, n):
            p, q = scipy.misc.pade(e_exp[:m + n + 1], n)
            return p.c[-1::-1], q.c[-1::-1]
        def print_result(p, q):
            optprint('num =', p)
            optprint('den =', q)
        def test_result(p, q, e_exp):
            m = len(p) - 1
            n = len(q) - 1
            # test against scipy
            p0, q0 = scipy_pade(m, n)
            assert numpy.allclose(gvar.mean(p), p0)
            assert numpy.allclose(gvar.mean(q), q0)
            # test that errors correlate with input coefficients
            num = gvar.powerseries.PowerSeries(p, order=m + n)
            den = gvar.powerseries.PowerSeries(q, order=m + n)
            ratio = (num/den).c / e_exp[:m + n + 1]
            assert numpy.allclose(gvar.mean(ratio), 1.)
            assert numpy.allclose(gvar.sdev(ratio), 0.0)

        # 1% noise --- automatically reduces to (2,1)
        e_exp_noise = [x * gvar.gvar('1.00(1)') for x in e_exp]
        p, q = pade_gvar(e_exp_noise, 3, 2, rtol=None)
        print_result(p, q)
        test_result(p, q, e_exp_noise)

        # 30% noise --- automatically reduces to (1,1)
        e_exp_noise = [x * gvar.gvar('1.0(3)') for x in e_exp]
        p, q = pade_gvar(e_exp_noise, 3, 2, rtol=None)
        print_result(p, q)
        test_result(p, q, e_exp_noise)

    def test_amu(self):
        " a_mu(vpol) "
        optprint('\n=========== Test a_mu')
        def no_vacpol(q2):
            return 0.25 / ALPHA ** 2

        # coefficient of alpha/pi
        amu = a_mu(no_vacpol)
        optprint('coef of alpha/pi = {}   error = {:7.2}'.format(
                amu, abs(amu-0.5) / 0.5
                ))
        assert numpy.allclose(amu, 0.5)

        # R. Karplus and N.M. Kroll result from Phys Rev 77 (#4), 536 (1950):
        # (alpha/pi)**2 * (3 + 11/36 - pi**2 / 3.)
        amu = a_mu(vacpol.fermion(m=Mmu))
        exact = (ALPHA/numpy.pi) ** 2 * ( 3 + 11./36 - numpy.pi**2 / 3.)
        optprint('a_mu(m=mu) = {}    error = {:7.2}'.format(
                amu, abs(amu/exact - 1.)
                ))
        assert numpy.allclose(amu/exact, 1.)

        # H. Suura and E.H. Wichmann in Phys Rev 105, 1930 (1950):
        # (alpha/pi)**2 ( log(mmu/me)/3 - 25/36 + O(me/mu))
        ratio = 1e5
        amu = a_mu(vacpol.fermion(Mmu/ratio))
        exact = (ALPHA/numpy.pi) ** 2 * ( numpy.log(ratio)/3. - 25./36.)
        assert numpy.allclose(amu/exact, 1., rtol=3/ratio)

    def test_noise(self):
        " a_mu(vpol) with noisy fermion loop "
        optprint('\n=========== Test noise (fermion loop)')
        def print_result(tag, amu, exact, pihat):
            line = '{:11} {:<13} {:15} {:15} {:15}'
            line = line.format(
                tag,
                amu if isinstance(amu, gvar.GVar) else '{:.8}'.format(amu),
                '  error = {:7.2}'.format(abs(gvar.mean(amu)/exact - 1.)),
                '  order = {}'.format(pihat.order),
                '  bad poles = {}'.format(pihat.badpoles())
                )
            optprint(line)

        # test at mK
        pihat_exact = vacpol.fermion(m=0.4937)
        exact = a_mu(pihat_exact)
        pihat = vacpol(pihat_exact.pseries['taylor'].c, (9,9))
        amu = a_mu(pihat)
        print_result('1loop(mK):', amu, exact, pihat)
        assert numpy.allclose(amu/exact, 1., rtol=1e-5)

        # mK with noise
        tayl = [
            ci * gvar.gvar('1.00(1)')
            for ci in pihat_exact.pseries['taylor'].c
            ]
        pihat = vacpol(tayl, (2,2))
        amu = a_mu(pihat)
        print_result('1loop(mK):', amu, exact, pihat)
        assert numpy.allclose(amu.mean/exact, 1., rtol=1e-2)

        # test at mpi
        pihat_exact = vacpol.fermion(m=MPI)
        exact = a_mu(pihat_exact)
        pihat = vacpol(pihat_exact.pseries['taylor'].c, (9,9))
        amu = a_mu(pihat)
        print_result('1loop(mpi):', amu, exact, pihat)
        assert numpy.allclose(amu/exact, 1., rtol=1e-4)

        # mpi with noise
        tayl = [
            ci * gvar.gvar('1.00(1)')
            for ci in pihat_exact.pseries['taylor'].c
            ]
        pihat = vacpol(tayl, (2,2))
        amu = a_mu(pihat)
        print_result('1loop(mpi):', amu, exact, pihat)
        assert numpy.allclose(amu.mean/exact, 1., rtol=1e-2)

    def test_scalar(self):
        " vacpole.scalar(mpi) "
        optprint('\n=========== Test scalar loop')
        for mpi, amu_vegas in [(MPI, '7.076903(1)e-9'), (MK, '6.631148(1)e-10')]:
            amu = a_mu(vacpol.scalar(mpi)) # a_mu_pipi(mpi)
            amu_vegas = gvar.gvar(amu_vegas)
            diff = gvar.fabs(amu - amu_vegas)
            assert diff.mean < 5 * diff.sdev
            optprint('1-loop({}) = {!s}   error = {}'.format(mpi, amu, diff))

    def test_exact_vs_pade(self):
        " a_mu from pade vs from function"
        optprint('\n=========== Test exact vs pade')
        m = MPI
        for n in [4, 5, 6, 7]:
            for f in [
                ('scalar', vacpol.scalar),
                ('fermion', vacpol.fermion),
                ('vector', vacpol.vector)
                ]:
                amu_exact = a_mu(f[1](m, n=n))
                vpol = f[1](m, n=n, use_pade=True)
                amu_pade = a_mu(vpol)
                optprint('{:>7}:  order = {}   pade/exact = {}'.format(
                    f[0], vpol.order, amu_pade / amu_exact
                    ))
                assert abs(amu_pade / amu_exact - 1.) < 0.01
            optprint(5 * '-')

if __name__ == '__main__':
    unittest.main()


"""
test_g2tools.py
"""
# Copyright (c) 2016-17 G. Peter Lepage.
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
        p0 = p0.c[::-1]
        q0 = q0.c[::-1]
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

    def test_pade_svd_consistency(self):
        " pade_svd self consistency "
        # high-order taylor series
        x = gv.powerseries.PowerSeries([0,1], order=20)
        f = np.exp(x).c
        # verify that reduced-order Pades are exact Pades
        m,n = 7,7
        for rtol in [1, 0.1, 0.01, 0.001]:
            a, b = pade_svd(f, m, n, rtol=rtol)
            mm = len(a) - 1
            nn = len(b) - 1
            if (m,n) != (mm,nn):
                aa, bb = pade_svd(f, mm, nn)
                self.assertTrue(np.allclose(aa, a))
                self.assertTrue(np.allclose(bb, b))

    def test_pade_gvar(self):
        " pade_gvar(tayl, m, n) "
        optprint('\n=========== Test pade_gvar')
        from scipy.misc import pade
        e_exp = [1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0, 1.0/720.]
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
            try:
                assert numpy.allclose(gvar.mean(p), p0)
            except:
                print (m,n, p0, p, q0, q)
            assert numpy.allclose(gvar.mean(q), q0)
            # test that errors correlate with input coefficients
            num = gvar.powerseries.PowerSeries(p, order=m + n)
            den = gvar.powerseries.PowerSeries(q, order=m + n)
            ratio = (num/den).c / e_exp[:m + n + 1]
            assert numpy.allclose(gvar.mean(ratio), 1.)
            assert numpy.allclose(gvar.sdev(ratio), 0.0)

        # print('scipy', scipy_pade(1,1), pade_svd(e_exp, 3,2, rtol=0.01))
        # 1% noise --- automatically reduces to (2,1)
        e_exp_noise = [x * gvar.gvar('1.0(1)') for x in e_exp]
        p, q = pade_gvar(e_exp_noise, 3, 2)
        print_result(p, q)
        self.assertEqual(len(p), 3)
        self.assertEqual(len(q), 2)
        test_result(p, q, e_exp_noise)

        # 30% noise --- automatically reduces to (1,1)
        e_exp_noise = [x * gvar.gvar('1.0(3)') for x in e_exp]
        p, q = pade_gvar(e_exp_noise, 3, 2)
        self.assertEqual(len(p), 2)
        self.assertEqual(len(q), 2)
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
        pihat = vacpol(pihat_exact.taylor(), (9,9))
        amu = a_mu(pihat)
        print_result('1loop(mK):', amu, exact, pihat)
        assert numpy.allclose(amu/exact, 1., rtol=1e-5)

        # mK with noise
        tayl = [
            ci * gvar.gvar('1.00(1)')
            for ci in pihat_exact.taylor()
            ]
        pihat = vacpol(tayl, (2,2))
        amu = a_mu(pihat)
        print_result('1loop(mK):', amu, exact, pihat)
        assert numpy.allclose(amu.mean/exact, 1., rtol=1e-2)

        # test at mpi
        pihat_exact = vacpol.fermion(m=MPI)
        exact = a_mu(pihat_exact)
        pihat = vacpol(pihat_exact.taylor(), (9,9))
        amu = a_mu(pihat)
        print_result('1loop(mpi):', amu, exact, pihat)
        assert numpy.allclose(amu/exact, 1., rtol=1e-4)

        # mpi with noise
        tayl = [
            ci * gvar.gvar('1.00(1)')
            for ci in pihat_exact.taylor()
            ]
        pihat = vacpol(tayl, (2,2), warn=True)
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

    def test_exact_vs_fourier(self):
        " a_mu from pade vs from function"
        optprint('\n=========== Test exact vs fourier')
        # fake data --- N=3 states
        N = 3
        ainv = 2.5
        Z = 1.5
        # the following are in lattice units, simulating lattice output
        m = np.array([0.5, 1.0, 1.5])[:N, None]
        t = np.arange(100)[None,:]
        G = np.sum(m / 4 * np.exp(-t*m), axis=0) / Z**2

        # fourier analysis
        fvpol = fourier_vacpol(G, ainv=ainv, Z=Z, periodic=False)
        a_mu_fourier = a_mu(fvpol, qmax=1000.)
        optprint('a_mu from fourier: {}'.format(a_mu_fourier))

        # exact result for 1, 2, and 3 states
        for n in range(1, N+1):
            a_mu_exact = np.sum(
                [a_mu(vacpol.vector(mi*ainv)) * ainv**2 for mi in m[:n]]
                )
            optprint('a_mu from {} states: {}'.format(n, a_mu_exact))
        self.assertLess(abs(1 - a_mu_fourier/a_mu_exact), 1e-4)

    def test_exact_vs_fourier_periodic(self):
        " a_mu from pade vs from function"
        optprint('\n=========== Test exact vs fourier')
        # loop over len(G) = even and odd
        for start in [-2, -1]:
            # fake data --- N=3 states
            N = 3
            ainv = 2.5
            Z = 1.5
            # the following are in lattice units, simulating lattice output
            m = np.array([0.5, 1.0, 1.5])[:N, None]
            t = np.arange(100)
            t = np.concatenate((t, t[start:0:-1]))
            G = np.sum(m / 4 * np.exp(-t*m), axis=0) / Z**2

            # fourier analysis
            fvpol = fourier_vacpol(G, ainv=ainv, Z=Z, periodic=True)
            a_mu_fourier = a_mu(fvpol, qmax=1000.)
            optprint('a_mu from fourier: {}'.format(a_mu_fourier))

            # exact result for 1, 2, and 3 states
            for n in range(1, N+1):
                a_mu_exact = np.sum(
                    [a_mu(vacpol.vector(mi*ainv)) * ainv**2 for mi in m[:n]]
                    )
                optprint('a_mu from {} states: {}'.format(n, a_mu_exact))
            self.assertLess(abs(1 - a_mu_fourier/a_mu_exact), 1e-4)

    def test_exact_vs_vacpol_FT(self):
        " a_mu from fourier_vacpol(vacpol.FT) "
        optprint('\n=========== Test exact vs fourier')
        # fake data --- N=3 states
        N = 3
        ainv = 2.5
        Z = 1.5
        # the following are in lattice units, simulating lattice output
        m = np.array([0.5, 1.0, 1.5])[:N, None]
        t = np.arange(100)[None,:]
        G = np.sum(m / 4 * np.exp(-t*m), axis=0) / Z**2
        t = t.reshape(-1)
        m = m.reshape(-1)

        # fourier analysis
        vpol = vacpol(moments(G, ainv=ainv, Z=Z, periodic=False), order=(3,3))
        self.assertLess(np.fabs(vpol.E[-1] - m[0] * ainv) / m[0]*ainv, 1e-6)
        self.assertLess(np.fabs(vpol.ampl[-1] - m[0] * ainv**3/4) / (m[0]*ainv**3/4), 1e-6)
        fvpol = fourier_vacpol(vpol.FT(t, ainv=ainv), ainv=ainv, periodic=False)
        a_mu_fmom = a_mu(fvpol, qmax=1000.)
        optprint('a_mu from FT of moments: {}'.format(a_mu_fmom))

        # exact result for 1, 2, and 3 states
        for n in range(1, N+1):
            a_mu_exact = np.sum(
                [a_mu(vacpol.vector(mi*ainv)) * ainv**2 for mi in m[:n]]
                )
            optprint('a_mu from {} states: {}'.format(n, a_mu_exact))
        self.assertLess(abs(1 - a_mu_fmom/a_mu_exact), 1e-4)

    def test_vacpol_poles(self):
        " vacpol.poles and vacpol.residues "
        m1 = gv.gvar('1.0(1)')
        f1 = gv.gvar('0.25(1)')
        vpol1 = vacpol.vector(m1, f=f1)
        m2 = gv.gvar('2.0(1)')
        f2 = gv.gvar('0.5(1)')
        vpol2 = vacpol.vector(m2, f=f2)
        # add two vectors together and check poles, residues
        vpol = vacpol(vpol1.taylor() + vpol2.taylor(), order=(2,2))
        self.assertEqual(gv.fabs(vpol.poles[0] + m2**2).fmt(5), '0.00000(0)')
        # print(gv.fabs(vpol.residues[0] + f2**2/2).fmt(5))
        self.assertEqual(gv.fabs(vpol.residues[0] + f2**2/2).fmt(5), '0.00000(0)')
        self.assertEqual(gv.fabs(vpol.poles[1] + m1**2).fmt(5), '0.00000(0)')
        self.assertEqual(gv.fabs(vpol.residues[1] + f1**2/2).fmt(5), '0.00000(0)')

    def test_warn_exception(self):
        " vacpol(warn=True) "
        # vacpol.scalar(MPI).taylor()
        tayl = np.array([  1.08361463e-02,  -3.97340348e-02,   2.26639708e-01,
            -1.58653778e+00,   1.25300529e+01,  -1.07205606e+02,
             9.71193290e+02,  -9.18408638e+03,   8.98033421e+04,
            -9.01971926e+05])
        tayl = tayl * gv.gvar(len(tayl) * ['1(1)'])
        tayl += np.array([  5.12534367e-05,  -2.70757996e-04,   5.49464167e-04,
            -2.69828134e-02,   9.43691955e-02,  -1.64530731e+00,
             4.97938388e-02,  -1.10418131e+01,  -7.24696697e+02,
             2.59030047e+04])
        with self.assertRaises(ValueError):
            vpol = vacpol(tayl, warn=True, qth=2*MPI, order=(3,3), rtol=1e-14)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # (1,1) Pade is exact for vector --- should reduce (3,3) to (1,1)
            m = gv.gvar('1.0(1)')
            f = gv.gvar('0.25(1)')
            vpol = vacpol(vacpol.vector(m, f, n=10).taylor(), order=(3,3), warn=True)
            self.assertTrue(w)
            self.assertEqual(vpol.order, (1,1))

if __name__ == '__main__':
    unittest.main()


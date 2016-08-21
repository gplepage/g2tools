""" Toolbox for muon g-2 analyses in lattice QCD.

This module contains a small number of tools useful for analyzing
contributions to the muon's magnetic moment from (lattice) QCD vacuum
polarization. The general technique is described in  arXiv:1403.1778.
The functions or classes include:

    moments(G)      --  compute moments of jj correlator G.
    mom2taylor(mom) --  convert moments in mom into Taylor coefficients
                        for q2-expansion.
    taylor2mom(tayl)--  convert Taylor expansion coefficients tayl into moments.
    vacpol(mom)     --  create a Pade approximant for the subtracted
                        vacuum polarization (Pi-hat) from the jj correlator
                        whose moments (or Taylor coefs) are in mom.
    a_mu(pihat, Q)  --  compute the contribution to the muon's g-2
                        anomaly from function pihat (usually built by vacpol).
    pade_gvar(f,m,n)--  general-purpose code for determining Pade approximants
                        to a power series whose coefficients are GVars.
    pade_svd(f,m,n) --  general-purpose code for determining Pade approximants
                        for a power series whose coefficients are floats.
                        Uses SVD regularization to stabilize results when
                        the input data are noisy.

A typical code sequence might look something like::

    mom = moments(G, ainv=ainv, periodic=True, nlist=[4, 6, 8, 10])
    vp = vacpol(mom, order=(2,2))
    amu = a_mu(vp, Q=Q)

where Q is the effective charge for the vacuum polarization (in units of
the proton charge). The first statement calculates the moments from G[t],
for n=4, 6, 8 and 10. The second statement constructs a vacuum polarization
function using the moments to create a (2,2) Pade approximant to Pi-hat.
The last statement calculates the contribution to the muon's anomaly
a_mu = (g-2)/2.

There are two equivalent representations used here to createa vacuum
polarization function. One is a dictionary G where G[n] = sum_t t**n G(t) is
the n-th moment (n=4,6,8...). The other is an array G of Taylor coefficients
G[i] where Pi-hat(q2) = q2 * sum_i q2**i G[i] (for i=0,1,2...). vacpol() can
take either type as an argument (dictionaries are moments, arrays are Taylor
coefficients).

Note that the following sample vacuum polarization functions are
included::

    vacpol.scalar(m)    --  scalar loop for scalar with mass m (eg, pion)
    vacpol.fermion(m)   --  fermion loop for fermion with mass m (eg, quark)
    vacpol.vector(m, f) --  tree-level for vector with mass m
                            and decay constant f (eg, rho)

This code requires the following Python modules: numpy, scipy, lsqfit, and gvar.
The last two are available on pypi and also at https://github.com/gplepage.

The code can be tested by running ``python g2tools.py``. Everything is
probably alright if there are no assertion errors. This code works with both
python2 (version > 2.7) and python3 (version > 3.4). Feel free to
contact g.p.lepage@cornell.edu if there are persistent problems.
"""

# Created by G. Peter Lepage (Cornell University) on 2014-09-13.
# Copyright (c) 2014-2016 G. Peter Lepage.
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

import gvar
import numpy
import lsqfit
import scipy.misc
import scipy.linalg
import scipy
import sys
import warnings
import collections

import math

__version__ = '1.0'

# constants
ALPHA = 1/137.035999074
Mmu = 0.1056583715
QMAX = 1e5
TOL = 1e-8
HMIN = 1e-12

def moments(G, Z=1., ainv=1., periodic=True, nlist=[4,6,8,10,12,14,16,18,20]):
    """ Compute t**n moments of correlator G.

    Compute ``sum_t t**n G(t)`` for ``n`` in ``nlist``, where both positive and
    negative ``t`` are included.

    Args:
        G: Array of correlator values ``G[t]`` for ``t=0,1...``.
        Z: Renormalization factor for current (moments multiplied by ``Z**2``).
            Defaul is 1.
        ainv: Inverse lattice spacing used to convert moments to
            physical units (n-th moment multiplied by ``1/ainv**(n-2)``).
            Default is 1.
        periodic: ``periodic=True`` implies ``G[-t] = G[t]`` (default);
            ``periodic=False`` implies no periodicity in array ``G[t]``
            (and results doubled to account for negative ``t``).
        nlist: List of moments to calculate. Default is
            ``nlist=[4,6,8...20]``.

    Returns:
        Dictionary ``Gmom`` where ``Gmom[n]`` is the ``n-th`` moment.
    """
    nt = len(G)
    if periodic:
        # sum over tmin <= t <= tmax
        tmax = math.floor(nt/2.)
        tmin = tmax - nt + 1.
        t = numpy.concatenate((numpy.arange(0., tmax + 1.), numpy.arange(tmin, 0.)))
        pfac = Z ** 2
    else:
        # sum over 0 <= t < nt and double, to account for t<0.
        t = numpy.arange(nt * 1.)
        pfac = 2. * Z ** 2
    Gmom = gvar.BufferDict()
    for n in nlist:
        Gmom[n] = numpy.sum( t ** n * G) * pfac / ainv ** (n-2)
    return Gmom

def a_mu(vacpol, Q=1, mmu=None, alpha=None, qmax=None, tol=None):
    """ Compute contribution to g-2 anomaly a_mu = (g-2)/2 from vacuum polarization.

    Args:
        vacpol: Function of ``q2`` for the subtracted vacuum polarization (Pi-hat).
            Here ``q2`` is space-like and so always positive.
            (See class :class:`vacpol`.)
        Q: Effective charge (in units of the proton's charge) --- for
            example, ``Q = 1./3. ``for s-quark loops (phi, etc) while
            ``Q = sqrt(5./9.)`` for u/d loops (rho, etc).
        mmu: Mass of the muon (default = Mmu).
        alpha: QED coupling (default = ALPHA).
        qmax: Maximum ``q`` included in integral (default is ``g2tools.QMAX = 1e5``).
        tol: Tolerance for integral over ``q2`` (default is ``g2tools.TOL = 1e-8``).

    Returns:
        Value of ``a_mu`` corresponding to ``Q**2 * vacpol``.
    """
    if mmu is None:
        mmu = Mmu
    if alpha is None:
        alpha = ALPHA
    if qmax is None:
        qmax = QMAX
    if tol is None:
        tol = TOL
    if hasattr(vacpol, 'badpoles') and vacpol.badpoles():
        raise RuntimeError('bad poles in vacpol: ' + str(vacpol.poles))
    fac = Q**2 * 4 * numpy.pi * alpha * (alpha / numpy.pi)
    def f_blum(q2):
        mmu2 = mmu**2
        z = - (q2 - numpy.sqrt(q2**2 + 4*q2*mmu2)) / 2 / q2 / mmu2
        return mmu2 * q2 * z**3 * (1 - q2 * z) / (1 + mmu2 * q2 * z**2)
    def deriv(theta, dummy):
        tantheta = numpy.tan(theta)
        q = tantheta * mmu
        jac = 2 * q * mmu * (tantheta ** 2 + 1.)
        q2 = q ** 2
        return jac * f_blum(q2) * vacpol(q2)
    odeint = gvar.ode.Integrator(deriv=deriv, tol=tol, hmin=HMIN)
    ans = odeint(0.0, interval=(1e-14, math.atan(qmax/mmu)))
    return ans * fac

def mom2taylor(mom):
    """ Convert moments in dictionary ``mom`` into Taylor series coefficients. """
    n = max(mom.keys()) // 2 - 1
    ans = n * [0.0]
    for i in sorted(mom):
        ans[i // 2 - 2] = (-1) ** (i/2) * mom[i] / math.factorial(i)
    return numpy.array(ans)

def taylor2mom(tayl):
    """ Convert Taylor coefficients in array ``tayl`` to moments. """
    ans = collections.OrderedDict()
    for i, ci in enumerate(tayl):
        n = 2 * i + 4
        ans[n] = (-1) ** (n/2) * math.factorial(n) * ci
    return ans

class vacpol(object):
    """ Subtracted vac. pol'n (``Pi-hat(q2)``) from correlator moments ``Gmon[n]``.

    The vacuum polarization function is a Pade approximant to the Taylor
    series corresponding to the moments ``g[n]``.  The code estimates the
    precision of the moments and sets the tolerance for the Pade determination
    accordingly. The order ``(m,n)`` of the Pade can be specified, but might
    be reduced by the code if the data are noisy.

    :class:`vacpol` objects are used primarily as functions (of *q*\ :sup:`2`)
    but also have several attributes. Attribute ``pseries`` is a dictionary
    containing various powerseries (see ``gvar.powerseries``) describing the
    function: the vacuum polarization function is ``q2`` times a Pade
    approximant  with a numerator given by ``pseries['num']`` and a
    denominator  given by ``pseries['den']``. The Taylor series for this
    function  is given by ``q2`` times ``pseries['taylor']``.

    :class:`vacpol` objects also have a method :func:`vacpol.badpoles` that
    tests the poles in the denomator of the Pade. ``badpoles(qth)`` returns
    ``False`` if any of the poles is complex or if any are located above
    ``-(qth ** 2)``. ``qth`` should be set equal to the threshold energy for
    the correlator. If it is unset, ``qth=0`` is used.

    :class:`vacpol` has several static methods for creating specialized
    examples of vacuum polarizations (e.g., for testing):

    * ``vacpol.fermion(m)`` -- 1-loop fermion (mass ``m``) contribution;

    * ``vacpol.scalar(m)`` -- 1-loop scalar (mass ``m``) contribution;

    * ``vacpol.vector(m, f)`` -- tree-level contribution from vector
        with mass ``m`` and decay constant ``f``.

    Args:
        g: Dictionary containing moments where ``g[n] = sum_t t**n G(t)``,
            or array containing Taylor coefficients where
            ``Pi-hat(q2) = q2 * sum_j q2**j * g[j]``.
        order: Tuple ``(m,n)`` specifying the order of the Pade
            approximant used to approximate the function. The order may
            be reduced (automatically) if the data are too noisy.
            If the order is not specified, it is set automatically
            according to the number of entries in ``G``.
        scale: Scale factor used to rescale ``q2`` so that
            the Taylor coefficients are more uniform in size. This is
            normally set automatically (from the first two moments),
            but the automatic value is overridden if ``scale`` is set.
        rtol: Relative tolerance assumed when determining the
            Pade approximant. This is normally set automatically
            (from the standard deviations of the moments), but the
            automatic value is overridden if ``rtol`` is specified.
    """
    def __init__(self, g, order=None, scale=None, rtol=None):
        f = mom2taylor(g) if hasattr(g, 'keys') else g
        if order is None:
            n = len(g)
            order = (n - n // 2, n // 2)
        m, n = order
        if scale is None:
            if len(f) > 1 and f[0] != 0 and f[1] != 0:
                scale = abs(gvar.mean(f[1] / f[0]))
            else:
                scale = 1.0
        self.scale = scale
        if scale != 1:
            fscaled = self.rescale(f, scale)
        else:
            fscaled = f
        p, q = pade_gvar(fscaled, m - 1, n, rtol=rtol)
        if scale != 1:
            p = self.rescale(p, 1./scale)
            q = self.rescale(q, 1./scale)
        self.fit = pade_gvar.fit
        self.rtol = pade_gvar.rtol
        self.pseries = dict(
            taylor=gvar.powerseries.PowerSeries(f),
            scaled_taylor=gvar.powerseries.PowerSeries(fscaled),
            num=gvar.powerseries.PowerSeries(p),
            den=gvar.powerseries.PowerSeries(q),
            )
        self.order = (len(p), len(q) - 1)
        self.poles = numpy.polynomial.polynomial.polyroots(gvar.mean(self.pseries['den'].c))

    @staticmethod
    def rescale(c, scale):
        " Rescale coefficients ``c``. "
        return c / (scale ** numpy.arange(len(c)))

    def __call__(self, q2):
        if hasattr(self, 'vacpol_fcn'):
            return self.vacpol_fcn(q2)
        return q2 * self.pseries['num'](q2) / self.pseries['den'](q2)

    def badpoles(self, qth=0):
        " False if any pole is complex or above threshold. "
        return numpy.any(numpy.iscomplex(self.poles)) or numpy.any(self.poles > -(qth ** 2))

    @staticmethod
    def fermion(m, n=19, use_pade=False):
        """ 1-loop subt. vac. pol'n from a fermion with mass m (and charge=1). """
        # m=1 taylor coefficients
        taylor_coef = numpy.array([
            0.0016886863940389628574, -0.00018093068507560316329,
            0.000026804545937126394562, -4.5689566938283627094e-6,
            8.4349969732215926943e-7, -1.6401383003486430239e-7,
            3.3078419502829775272e-8, -6.8550671995995915859e-9,
            1.4508078729311304944e-9, -3.1223908569604764988e-10,
            6.8124891424592214518e-11, -1.5033733755735627587e-11,
            3.3496913408005111864e-12, -7.5252282425817935639e-13,
            1.7026779053922442003e-13, -3.8766327310269845632e-14,
            8.8749620551333191909e-15, -2.0417469685242607540e-15,
            4.7177721994142097012e-16
            ], object)
        n = min(len(taylor_coef), n)
        order = (n - n // 2, n // 2)
        vpol = vacpol(
            taylor_coef[:n] * (1. / m ** 2) ** (1 + numpy.arange(n)),
            order=order,
            )
        def vacpol_fcn(q2):
            """ from Lifshitz and Pitaevskii vol 2

            Need to treat xi near 1 differently because of
            roundoff error. Use the power series (in q2/m2)
            in that region, instead of the exact formula.
            The boundary between the two regions is xi_bdy.
            Varying that (from 0.2 to 0.95) confirms the
            consistency of the power series with the formula.
            """
            xi_bdy = 0.8    # 1-xi is q/m roughly
            q2m2 = q2 / m ** 2
            xi = 2. / (2. + q2m2 + numpy.sqrt(q2m2 ** 2 + 4 * q2m2))
            if xi < xi_bdy:
                # use formula except for q2 near zero
                ans = - 1./12./numpy.pi**2 * xi / (1 - xi) ** 2 * (
                    - 22. / 3.
                    + 5. / 3. * (xi + 1. / xi)
                    + (xi + 1. / xi - 4.) * (1 + xi)
                    / (1 - xi) * numpy.log(xi)
                    )
            else:
                # use taylor expansion for q2 close to 0
                ans = 0.0
                for n, cn in enumerate(taylor_coef):
                    ans += cn * q2m2 ** (n+1)
            return ans
        if not use_pade:
            vpol.vacpol_fcn = vacpol_fcn
        return vpol

    @staticmethod
    def scalar(m, n=10, use_pade=False):
        """ 1-loop subt. vac. pol'n from a scalar with mass m (and charge=1). """
        j = numpy.arange(n) + 1.
        fact = scipy.misc.factorial
        taylor_coef = (1/ 8. / numpy.pi**2) * (-1) ** (j+1) / m ** (2 * j) * (
            fact(j + 1) * fact(j - 1) * 1. / fact(2 * j + 3)
            )
        order = (n - n // 2, n // 2)
        vpol = vacpol(taylor_coef, order=order)
        def h(q2):
            y = 4 * m**2 / q2
            sqrt1_y = (1 + y) ** 0.5
            G = 0.5 / sqrt1_y * numpy.log((1 + sqrt1_y) ** 2 / y)
            return 2./3. + 2 * (1+y) - 2 * (1+y)**2 * G
        def vacpol_fcn(q2):
            if q2 < m ** 2 * 1e-3:
                ans = vpol.pseries['taylor'](q2) * q2
            else:
                ans = -h(q2) / 48. / numpy.pi ** 2
            return ans
        if not use_pade:
            vpol.vacpol_fcn = vacpol_fcn
        return vpol

    @staticmethod
    def vector(m, f=1, n=10, use_pade=False):
        """ Vac. pol'n due to a vector with mass ``m`` and decay const. ``f``.

        The decay constant is defined such that the vacuum polarization
        function is ``Pi-hat = q2 * f**2/2/m**2 / (q2 + m**2)``.
        """
        n = 20
        j = numpy.arange(n) + 1.
        taylor_coef = f ** 2 / 2. / m ** (2 * j + 2) * (-1) ** (j + 1)
        order = (n - n // 2, n // 2)
        vpol = vacpol(taylor_coef, order=order)
        def vacpol_fcn(q2):
            return q2 * f ** 2 / (q2 + m ** 2) / 2. /m**2
        if not use_pade:
            vpol.vacpol_fcn = vacpol_fcn
        return vpol

def pade_gvar(f, m, n, rtol=None):
    """ ``[m,n]`` Pade approximant to ``sum_i f[i] x**i`` for ``GVar``\s.

    The ``[m,n]`` Pade approximant to a series given by
    ``sum_i f[i] * x**i`` is the ratio of  polynomials of order ``m``
    (numerator) and ``n`` (denominator) whose  Taylor expansion agrees
    with that of the original series up to order ``m+n``.

    This code uses an *svd* algorithm (see :func:`pade_svd`) to deal with
    imprecision in the input data. It automatically reduces
    the order of the approximant if the extraction of Pade coefficients
    is too unstable given noise in the input data.

    :param f: Array ``f[i]`` of power series coefficients for ``i=0...n+m``.
    :param m: Maximum order of polynomial in numerator of Pade
        approximant (``m>=0``).
    :param n: Maximum order of polynomial in denominator of Pade
        approximant (``m>=0``).
    :param rtol: Relative accuracy of input coefficients. Overrides
        default estimate from the ``f[i]`` unless set equal to ``None``.
    :returns: Tuple of power series coefficients ``(p, q)`` such that
        ``sum_i p[i] x**i`` is the numerator of the approximant,
        and ``sum_i q[i] x**i`` is the denominator. ``q[0]`` is
        normalized to 1.
    """
    # check inputs
    if not numpy.any([isinstance(fi, gvar.GVar) for fi in f]):
        pade_gvar.fit = None
        pade_gvar.rtol = rtol
        return pade_svd(f, m, n, rtol=1e-14 if rtol is None else rtol)
    c = numpy.array(f[:n + m + 1])
    if len(f) < (m + n + 1):
        raise ValueError(
            'not enough f[i]s -- need {} have {}'.format(n + m + 1, len(f))
            )

    # compute tolerance if not specifiec
    if rtol is None:
        rtol = gvar.sdev(numpy.sum(c))
        if rtol == 0:
            rtol = 1e-14
        else:
            rtol /= numpy.sum(numpy.abs(gvar.mean(c)))

    # find approximate means
    p, q = pade_svd(gvar.mean(c), m, n, rtol=rtol)
    p0 = dict(num=p, den=q[1:])
    m = len(p) - 1
    n = len(q) - 1

    # fit to insert errors
    def fitfcn(p):
        order = m + n
        num = gvar.powerseries.PowerSeries(p['num'], order=order)
        den = gvar.powerseries.PowerSeries([1.] + p['den'].tolist(), order=order)
        ratio = num / den
        return ratio.c
    fit = lsqfit.nonlinear_fit(
        data=c[:n + m + 1],
        fcn=fitfcn,
        p0=p0,
        reltol=1e-10,
        debug=True
        )
    if fit.chi2 > 1.:
        warnings.warn('bad fit: chi2 = {}'.format(fit.chi2))

    # save intermediate results in case they are needed later
    pade_gvar.fit = fit
    pade_gvar.rtol = rtol
    pade_gvar.p0 = p0
    return fit.p['num'], numpy.array([1.] + fit.p['den'].tolist())

def pade_svd(f, m, n, rtol=1e-14):
    """ ``[m,n]`` Pade approximant to ``sum_i f[i] x**i``.

    The ``[m,n]`` Pade approximant to a series given by
    ``sum_i f[i] * x**i`` is the ratio of  polynomials of order ``m``
    (numerator) and ``n`` (denominator) whose  Taylor expansion agrees
    with that of the original series up to order ``m+n``.

    This code is adapted from P. Gonnet,  S. Guttel, L. N. Trefethen, SIAM
    Review Vol 55, No. 1, 101 (2013). It uses an *svd* algorithm to deal with
    imprecision in the input data,  here specified by the relative tolerance
    ``rtol`` for the  input coefficients ``f[i]``. It automatically reduces
    the order of the approximant if the extraction of Pade coefficients
    is too unstable given tolerance ``rtol``.

    :param f: Array ``f[i]`` of power series coefficients for ``i=0...n+m``.
    :param m: Maximum order of polynomial in numerator of Pade
        approximant (``m>=0``).
    :param n: Maximum order of polynomial in denominator of Pade
        approximant (``m>=0``).
    :param rtol: Relative accuracy of input coefficients.
    :returns: Tuple of power series coefficients ``(p, q)`` such that
        ``sum_i p[i] x**i`` is the numerator of the approximant,
        and ``sum_i q[i] x**i`` is the denominator. ``q[0]`` is
        normalized to 1.
    """
    linalg = scipy.linalg
    c = numpy.array(f[:n + m + 1], float)
    if len(f) < (m + n + 1):
        raise ValueError(
            'not enough f[i]s -- need {} have {}'.format(n + m + 1, len(f))
            )
    ts = rtol * linalg.norm(c)
    if linalg.norm(c[:m + 1]) <= rtol * linalg.norm(c):
        a = numpy.array([0.])
        b = numpy.array([1.])
        mu = None
        nu = 0
    else:
        row = numpy.zeros(n+1)
        row[0] = c[0]
        col = c
        while True:
            if n == 0:
                a = c[:m+1]
                b = numpy.array([1.])
                return a, b
            Z = linalg.toeplitz(col[:m + n + 1], row[:n + 1])
            C = Z[m + 1:, :]
            rho = numpy.sum(linalg.svdvals(C) > ts)
            if rho == n:
                break
            m -= n - rho
            n = rho
        if n > 0:
            # use svd to get solution b, but only to normalize C
            # then use QR decomposition to get final b
            U, S, V = linalg.svd(C, full_matrices=True)
            b = V.transpose()[:, -1]
            D = numpy.diag(numpy.abs(b) + numpy.sqrt(sys.float_info.epsilon))
            Q,R = linalg.qr(C.dot(D).transpose())
            b = D.dot(Q)[:,-1]
            b = b / linalg.norm(b)
            a = Z[:m + 1, :n + 1].dot(b)
            lam = numpy.where(abs(b) > rtol)[0][0]
            b = b[lam:]
            a = a[lam:]
            b = b[:numpy.where(abs(b) > rtol)[0][-1] + 1]
        try:
            a = a[:numpy.where(abs(a) > ts)[0][-1] + 1]
        except IndexError:
            a = a[:1]
        a = a/b[0]
        b = b/b[0]
    return a, b

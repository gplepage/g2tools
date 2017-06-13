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
    fourier_vacpol(G)-- create subtracted vacuum polarization (``PI-hat``) by
                        Fourier transforming *jj* correlator ``G(t)``.
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
"""

# Created by G. Peter Lepage (Cornell University) on 2014-09-13.
# Copyright (c) 2014-2017 G. Peter Lepage.
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

__version__ = '1.3'

USE_SCIPY_PADE = False

# constants
ALPHA = 1/137.035999074
Mmu = 0.1056583715
QMIN = 1e-15
QMAX = 1e5
TOL = 1e-8
HMIN = 1e-12

def moments(G, Z=1., ainv=1., periodic=True, tmin=None, nlist=[4,6,8,10,12,14,16,18,20]):
    """ Compute t**n moments of correlator G.

    Compute ``sum_t t**n G(t)`` for ``n`` in ``nlist``, where both positive and
    negative ``t`` are included.

    Args:
        G: Array of correlator values ``G[t]`` for ``t=0,1...`` (in
            lattice units).
        Z: Renormalization factor for current (moments multiplied by ``Z**2``).
            Defaul is 1.
        ainv: Inverse lattice spacing used to convert moments to
            physical units (n-th moment multiplied by ``1/ainv**(n-2)``).
            Default is 1.
        periodic: ``periodic=True`` implies ``G[-t] = G[t]`` (default);
            ``periodic=False`` implies no periodicity in array ``G[t]``
            (and results doubled to account for negative ``t``).
        tmin: minimum ``t`` value included in moments; ignored
            if ``None`` (default).
        nlist: List of moments to calculate. Default is
            ``nlist=[4,6,8...20]``.

    Returns:
        Dictionary ``Gmom`` where ``Gmom[n]`` is the ``n-th`` moment.
    """
    nt = len(G)
    if periodic:
        # sum over t1 <= t <= t2
        t2 = math.floor(nt/2.)
        t1 = t2 - nt + 1.
        t = numpy.concatenate((numpy.arange(0., t2 + 1.), numpy.arange(t1, 0.)))
        pfac = Z ** 2
    else:
        # sum over 0 <= t < nt and double, to account for t<0.
        t = numpy.arange(nt * 1.)
        pfac = 2. * Z ** 2
    Gmom = gvar.BufferDict()
    if tmin is not None:
        idx = numpy.fabs(t) >= tmin
        t = t[idx]
        G = G[idx]
    for n in nlist:
        Gmom[n] = numpy.sum( t ** n * G) * pfac / ainv ** (n-2)
    return Gmom

def a_mu(
    vacpol, Q=1, mmu=None, alpha=None, qmin=None, qmax=None, rescale=None,
    tol=None, exceptions=True
    ):
    """ Compute contribution to g-2 anomaly a_mu = (g-2)/2 from vacuum polarization.

    Args:
        vacpol: Function of ``q2`` for the subtracted vacuum polarization (Pi-hat).
            Here ``q2`` is space-like and so always positive.
            (See class :class:`vacpol`.)
        Q: Effective charge (in units of the proton's charge) --- for
            example, ``Q = 1./3. ``for s-quark loops (phi, etc) while
            ``Q = sqrt(5./9.)`` for u/d loops (rho, etc). (Default is 1.)
        mmu: Mass of the muon (default is ``g2tools.Mmu``).
        alpha: QED coupling (default is ``g2tools.ALPHA``).
        qmin: Maximum ``q`` included in integral (default is ``g2tools.QMIN = 1e-15``).
        qmax: Maximum ``q`` included in integral (default is ``g2tools.QMAX = 1e5``).
        rescale: Rescales momentum in vacuum poln: ``vacpol(q2 * rescale**2)``
            (default is 1).
        tol: Tolerance for integral over ``q2`` (default is ``g2tools.TOL = 1e-8``).
        exceptions: If ``True`` (default), an exception is raised if there
            are bad poles in the ``vacpol``. If ``False``, exceptions are
            suppressed.

    Returns:
        Value of ``a_mu`` corresponding to ``Q**2 * vacpol``.
    """
    if mmu is None:
        mmu = Mmu
    if alpha is None:
        alpha = ALPHA
    if qmin is None:
        qmin = QMIN
    if qmax is None:
        qmax = QMAX
    if tol is None:
        tol = TOL
    if rescale is None:
        r2 = 1.
    r2 = 1. if rescale is None else rescale * rescale
    if exceptions and hasattr(vacpol, 'badpoles') and vacpol.badpoles():
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
        return jac * f_blum(q2) * vacpol(q2*r2)
    odeint = gvar.ode.Integrator(deriv=deriv, tol=tol, hmin=HMIN)
    ans = odeint(0.0, interval=(math.atan(qmin/mmu), math.atan(qmax/mmu)))
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

class fourier_vacpol(object):
    """ Subtracted vac. pol'n (``Pi-hat(q2)``) from correlator ``G(t)``.

    The correlator is Fourier transformed to produce a function ``Pi_hat``
    of (Euclidean) *q*\ :sup:`2` suitable for use in :func:`g2tools.a_mu`.

    See Bernecker & Meyer, EPJA47 (2011) 148 , arXiv:1107.4388 for details
    on the Fouier transformation.

    Args:
        G (array): Current-current correlator in an array whose elements
            are ``[G(0),G(a),G(2*a),...,G(-2*a),G(-a)]`` if
            ``periodic=True`` or ``[G(0),G(a),...,G(T*a-1)]``
            otherwise. ``G`` is assumed to be in lattice units.
        Z: Renormalization factor for current (correlator multiplied
            by ``Z**2``). Defaul is 1.
        ainv: Inverse lattice spacing used to convert Fourier transform to
            physical units. Default is 1.
        periodic: ``periodic=True`` implies ``G[-t] = G[t]`` (default);
            ``periodic=False`` implies ``G[t]`` is not periodic and
            is specified for only non-negative ``t`` values
            (results are doubled to account for negative ``t``).
    """
    def __init__(self, G, Z=1., ainv=1., periodic=True):
        G = numpy.array(G)
        if periodic:
            nG = len(G)
            if nG % 2 == 0:
                self.G = G[:nG // 2 + 1]
                self.G[1:-1] += G[-1:nG // 2:-1]
                self.G[1:-1] *= 0.5
            else:
                self.G = G[:nG//2 + 1]
                self.G[1:] += G[-1:nG//2:-1]
                self.G[1:] *= 0.5
        else:
            self.G = G
        # In next expression should be ainv**3 except for factor of 1/ainv
        # associated with the t integral in __call__. So self.G_ainv is
        # G / ainv.
        self.G_ainv = self.G * Z**2 * ainv**2
        self.t = numpy.arange(len(self.G_ainv)) / ainv   # q is in phys units
        self.t2 = self.t ** 2

    def __call__(self, q2):
        # N.B. factor of 1/ainv needed for t-integral inncluded in self.G
        return numpy.sum(
            (self.t2 - 4 * (numpy.sin((q2 ** 0.5/2.) * self.t))**2 / q2) *
            self.G_ainv
            )

class vacpol(object):
    """ Subtracted vac. pol'n (``Pi-hat(q2)``) from correlator moments ``Gmon[n]``.

    The current-current correlator is ``q2 * Pi(q2)``, where
    ``Pi-hat(q2) = Pi(q2) - Pi(0)`` is the subtracted (i.e., renormalized)
    vacuum polariztion function.

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
    the correlator. If it is unset, ``qth=0`` is used. Lists of the
    poles and their residues (for ``Pi-hat(q2)``) are available in
    attributes ``pole`` and ``residue``, respectively.

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
        qth: Threshold for particle production: poles above ``-qth**2``
            are bad. Default is ``qth=0``.
        warnings: ``warnings=True`` causes a warning to be issued when the
            order has been reduced automatically. ``warnings=False`` (default)
            suppresses the warnings.
        exceptions: If ``True`` (default), an exception is raised if there
            are bad poles in the ``vacpol``. If ``False``, exceptions are
            suppressed.
    """
    def __init__(
        self, g, order=None, scale=None, rtol=None, qth=0,
        warn=True, exceptions=True
        ):
        f = mom2taylor(g) if hasattr(g, 'keys') else g
        self.qth = qth
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
        # self.fit = pade_gvar.fit
        # self.rtol = pade_gvar.rtol
        self.pseries = dict(
            taylor=gvar.powerseries.PowerSeries(f),
            scaled_taylor=gvar.powerseries.PowerSeries(fscaled),
            num=gvar.powerseries.PowerSeries(p),
            den=gvar.powerseries.PowerSeries(q),
            )
        self.order = (len(p), len(q) - 1)
        self.poles = numpy.polynomial.polynomial.polyroots(
            gvar.mean(self.pseries['den'].c)
            )
        if (exceptions and (
            numpy.any(numpy.iscomplex(self.poles)) or
            numpy.any(self.poles > -(self.qth ** 2))
            )):
            raise ValueError(
                'bad poles (try reducing the order): ' +
                str(self.poles)
                )
        if warn == True:
            if self.order != order:
                warnings.warn('reduced order: {}'.format(self.order))
            if self.badpoles():
                warnings.warn('bad poles: {}', str(self.poles))
        # Fourier parameters
        num = self.pseries['num']
        den = self.pseries['den']
        res = []
        pole = []
        # add errors to poles and residues
        def add_sdev(poly, p):
            if not isinstance(poly(p), gvar.GVar):
                return p
            nc = gvar.mean(poly.c) + gvar.gvar(len(poly.c) * ['0(0)'])
            np = p + gvar.gvar('0(0)')
            npoly =  gvar.powerseries.PowerSeries(nc)
            npoly_np = npoly(np)
            dpoly_dp = npoly_np.deriv(np)
            for (nci,ci) in zip(npoly.c, poly.c):
                if not isinstance(ci, gvar.GVar):
                    continue
                p = p - npoly_np.deriv(nci) * (ci - ci.mean) / dpoly_dp
            return p
        for p in self.poles:
            pole.append(add_sdev(den, p))
            res.append(num(pole[-1]) * pole[-1] / den.c[-1])
        for i,p in enumerate(pole):
            other_p = numpy.array(pole[:i] + pole[i+1:])
            res[i] /= numpy.prod(p-other_p)
        # residue is residue of pole in Pi-hat(q2) (not q2 * Pi-hat(q2))
        self.residues = numpy.array(res)
        self.poles = numpy.array(pole)
        self.E = (-self.poles) ** 0.5
        # ampl is amplitude for q2 * Pi-hat(q2) hence extra factor of pole.
        self.ampl = self.residues * self.poles/ 2. / self.E
        # calculate residual polynomial - not used for anything!
        # q2 = gvar.powerseries.PowerSeries(
        #     [0., 1.], order=num.order - den.order + 1
        #     )
        # self.direct = (
        #     q2 * num / den -
        #     numpy.sum( [r / (q2 - p) for r,p in zip(self.residues, self.poles)])
        #     )

    @staticmethod
    def rescale(c, scale):
        " Rescale coefficients ``c[j] -> c[j] / scale**j`` "
        return c / (scale ** numpy.arange(len(c)))

    def __call__(self, q2):
        if hasattr(self, 'vacpol_fcn'):
            return self.vacpol_fcn(q2)
        return q2 * self.pseries['num'](q2) / self.pseries['den'](q2)

    def taylor(self, n=None):
        """ Return Taylor coefficients for ``PI-hat(q2)/q2``.

        Args:
            n: Maximum number of coefficients returned. Returns
                all coefficents if ``None`` (default)/
        """
        return numpy.array(self.pseries['taylor'].c[:n])

    def badpoles(self, qth=None):
        """ True if any pole is complex or above threshold.

        Args:
            qth: Threshold for particle production: poles above ``-qth**2``
                are bad. (Default is ``qth=0``.)
        """
        if qth is None:
            qth = self.qth
        return numpy.any(numpy.iscomplex(self.poles)) or numpy.any(self.poles > -(qth ** 2))

    def FT(self, t, ainv=1.):
        """ Fourier transform of ``q2 * PI-hat(q2)``.

        The Pade approximant can be decomposed into a sum of poles (partial
        fractions), which give a sum of decaying exponentials when Fourier
        transformed back to t-space. The amplitudes and energies of these
        exponentials (for the transform of ``q2 * Pi-hat(q2)'') are stored in
        :class:`g2tools.vacpol` attributes ``E`` and ``ampl``, respectively.

        The decomposition into a sum of poles leaves a residual polynomial in
        ``q2`` (zeroth-order for ``(n,n)`` Pades). This is ignored in the
        Fourier transform since it typically affects the transform only for
        very small ``t``. These terms have a negligible effect (suppressed by
        ``a**2j`` on the Taylor coefficients ``Pi[j]`` of ``Pi-hat(q2)``
        (for j>=1).

        Optional parameter ``ainv`` can be used to convert the Fourier
        transform to lattice units (by multiplying it by ``1/ainv**3``) for
        comparison with simulation data. The times ``t`` are then
        assumed to be in lattice units.

        Args:
            t (number, array): Time in physical units unless ``ainv``
                is specified, in which case lattice units are assumed.
            ainv: Inverse lattice spacing. The Fourier transform is in lattice
                units if ``ainv`` is specified (assuming the original
                Taylor coefficients are in physical units).
        """
        # Need 1/ainv**3 to put G(t) into lattice units.
        if numpy.shape(t) == ():
            return numpy.sum(self.ampl * numpy.exp(-self.E*t)) / ainv ** 3
        else:
            tshape = numpy.shape(t)
            t = numpy.asarray(t).flatten() / ainv
            ans = numpy.sum(
                self.ampl[:, None] * numpy.exp(-self.E[:, None]*t[None, :]),
                axis=0
                )
            return ans.reshape(tshape)/ ainv ** 3

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
    def vector(m, f=1., n=10, use_pade=False):
        """ Vac. pol'n due to a vector with mass ``m`` and decay const. ``f``.

        The decay constant is defined such that the vacuum polarization
        function is ``Pi-hat = q2 * f**2/2/m**2 / (q2 + m**2)``.
        """
        j = numpy.arange(n) + 1.
        taylor_coef = f ** 2 / 2. / m ** (2 * j + 2) * (-1) ** (j + 1)
        order = [(0,0), (1,0), (1,1)][min(2, n)]
        vpol = vacpol(taylor_coef, order=order)
        def vacpol_fcn(q2):
            return q2 * f ** 2 / (q2 + m ** 2) / 2. /m**2
        if not use_pade:
            vpol.vacpol_fcn = vacpol_fcn
        return vpol

def pade_gvar(f, m, n, rtol='gavg'):
    """ ``[m,n]`` Pade approximant to ``sum_i f[i] x**i`` for ``GVar``\s.

    The ``[m,n]`` Pade approximant to a series given by
    ``sum_i f[i] * x**i`` is the ratio of  polynomials of order ``m``
    (numerator) and ``n`` (denominator) whose  Taylor expansion agrees
    with that of the original series up to order ``m+n``.

    This code uses an SVD algorithm (see :func:`pade_svd`) to deal with
    imprecision in the input data. It automatically reduces
    the order of the approximant if the extraction of Pade coefficients
    is too unstable given noise in the input data.

    Args:
        f: Array ``f[i]`` of power series coefficients for ``i=0...n+m``.
        m: Maximum order of polynomial in numerator of Pade
            approximant (``m>=0``).
        n: Maximum order of polynomial in denominator of Pade
            approximant (``m>=0``).
        rtol (float or str): If ``rtol`` is a string, it determines how the
            relative tolerance is determined from the relative
            uncertainties in the ``f[i]``. Set ``rtol`` equal to:
            ``'gavg'`` for the geometric mean (default); ``'avg'`` for
            the average; ``'min'`` for the minimum; or ``'max'`` for
            the maximum. Otherwise a number can be specified, in which case
            the uncertainties in ``f[i]`` are ignored.
    Returns:
        Tuple of power series coefficients ``(p, q)`` such that
        ``sum_i p[i] x**i`` is the numerator of the approximant,
        and ``sum_i q[i] x**i`` is the denominator. ``q[0]`` is
        normalized to 1.
    """
    # check inputs
    assert m >= 0 and n >= 0
    f = f[:n + m + 1]
    if len(f) < (m + n + 1):
        raise ValueError(
            'not enough f[i]s -- need {} have {}'.format(n + m + 1, len(f))
            )
    if not numpy.any([isinstance(fi, gvar.GVar) for fi in f]):
        pade_gvar.rtol = 1e-14
        return pade_svd(f, m, n)
    else:
        c = numpy.array(f)

    # compute tolerance if not specified
    if rtol in ['avg', 'min', 'max', 'gavg']:
        means = numpy.fabs(gvar.mean(c))
        sdevs = gvar.sdev(c)
        idx = means > 0.0
        if numpy.any(idx) and numpy.all(sdevs[idx] > 0):
            ratio = sdevs[idx] / means[idx]
            if rtol == 'gavg':
                # geometric mean
                rtol = numpy.exp(
                    numpy.average(numpy.log(ratio))
                    )
            elif rtol == 'avg':
                rtol = numpy.average(ratio)
            elif rtol == 'min':
                rtol = numpy.min(ratio)
            else:
                rtol = numpy.max(ratio)
        else:
            rtol = 1e-14
    elif rtol is not None:
        rtol = numpy.fabs(rtol)
    else:
        rtol = 1e-14
    pade_gvar.rtol = rtol

    # find Pade coefficientw
    p, q = pade_svd(gvar.mean(c), m, n, rtol=rtol)
    m = len(p) - 1
    n = len(q) - 1

    # add uncertainties
    p = p * gvar.gvar(len(p) * ['1(0)'])
    q = q[1:] * gvar.gvar(len(q[1:]) * ['1(0)'])
    num = gvar.powerseries.PowerSeries(p, order=m + n)
    den = gvar.powerseries.PowerSeries([1] + list(q), order=m + n)
    pade = numpy.concatenate((p,q))
    cc = (num / den).c
    M = numpy.empty((len(pade), len(pade)), float)
    for i in range(len(pade)):
        for j in range(len(pade)):
            M[i, j] = cc[i].deriv(pade[j])
    pade = pade + gvar.linalg.solve(M, (c - gvar.mean(c))[:len(pade)])
    return (
        numpy.array(pade[:m + 1]),
        numpy.array([gvar.gvar(1,0)] + list(pade[m + 1:]))
        )


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
    :param rtol: Relative accuracy of input coefficients. (Default is 1e-14.)
    :returns: Tuple of power series coefficients ``(p, q)`` such that
        ``sum_i p[i] x**i`` is the numerator of the approximant,
        and ``sum_i q[i] x**i`` is the denominator. ``q[0]`` is
        normalized to 1.
    """
    linalg = scipy.linalg
    mn_save = m,n
    c = numpy.array(f[:n + m + 1], float)
    if len(f) < (m + n + 1):
        raise ValueError(
            'not enough f[i]s -- need {} have {}'.format(n + m + 1, len(f))
            )
    if USE_SCIPY_PADE:
        p, q = scipy.misc.pade(c, n)
        return numpy.array(p.c[::-1]), numpy.array(q.c[::-1])
    ts = rtol * linalg.norm(c)
    if linalg.norm(c[:m + 1]) <= rtol * linalg.norm(c):
        # return power series through order m
        a = numpy.array(c[:1])
        b = numpy.array([1.])
    else:
        row = numpy.zeros(n+1)
        row[0] = c[0]
        col = c
        while True:
            if n == 0:
                # return the power series through order m
                a = c[:m + 1]
                b = numpy.array([1.])
                return a, b
            Z = linalg.toeplitz(col[:m + n + 1], row[:n + 1])
            C = Z[m + 1:, :]
            rho = numpy.sum(linalg.svdvals(C) > ts)
            if rho == n:
                break
            m -= n - rho
            n = rho
            if m < 0:
                m = 0
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
        idx = abs(a) > ts
        if not numpy.any(idx):
            a = a[:1]
        else:
            a = a[:numpy.where(idx)[0][-1] + 1]
        a = a / b[0]
        b = b / b[0]
    # Ending modified so that an approximant for non-zero rtol
    # is the same as the reduced-order approximant evaluated with
    # zero rtol; any approximant returned by the algorithm should
    # be an exact approximant to the input. Thus rtol determines the
    # order of the final approximant, but does not affect the values
    # of the approximant's coefficients. This is not true of the original.
    mfinal = len(a) - 1
    nfinal = len(b) - 1
    return (a,b) if (mfinal,nfinal) == mn_save else pade_svd(f, mfinal, nfinal)

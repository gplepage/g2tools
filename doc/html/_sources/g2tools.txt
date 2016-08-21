:mod:`g2tools` Module
=====================================================

.. module:: g2tools
   :synopsis: Least-Squares Fit to Correlators.

.. moduleauthor:: G.P. Lepage <g.p.lepage@cornell.edu>

.. |g2tools| replace:: :mod:`g2tools`

.. |GVar| replace:: :class:`gvar.GVar`

Moments
----------
The main tools for creating and manipulating moments are:

.. autofunction:: g2tools.moments

.. autofunction:: g2tools.mom2taylor

.. autofunction:: g2tools.taylor2mom

Subtracted Vacuum Polarization
-------------------------------
A subtracted vacuum polarization function (``Pi-hat``) is
represented by the following class:

.. autoclass:: g2tools.vacpol


Padé Approximants
-------------------

The following two functions are used for calculating Padé approximants from
the Taylor coefficients of an arbitrary function. The first
(:func:`g2tools.pade_svd`) implements an algorithm that uses *svd* cuts to
address instabilities caused  by uncertainties in the Taylor coefficients. The
second function (:func:`g2tools.pade_gvar`) is built on the first but allows
Taylor coefficients to have uncertainties (|GVar|\s). The statistical
uncertainties and correlations between different coefficients are propagated
through the analysis.

.. autofunction:: g2tools.pade_svd

.. autofunction:: g2tools.pade_gvar

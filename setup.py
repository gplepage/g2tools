from distutils.core import setup
import g2tools

setup(name='g2tools',
    version=g2tools.__version__,
    description='Utilities for muon g-2 analyses in lattice QCD.',
    author='G. Peter Lepage, Cornell University',
    author_email='g.p.lepage@cornell.edu',
    license='GPLv3+',
    py_modules=['g2tools'],
    requires=["lsqfit (>=9.1)", 'numpy (>=1.7)', 'gvar (>=7.3)', 'scipy'],
    install_requires=['lsqfit>=9.1', 'gvar>=7.3', 'numpy>=1.7', 'scipy'],
    platforms="Any",
    url="https://github.com/gplepage/g2tools.git",
    long_description="""\
    This module contains a small number of tools useful for analyzing
    contributions to the muon's magnetic moment from (lattice) QCD vacuum
    polarization. The general technique is described in  arXiv:1403.1778.
    The functions or classes include:

        moments(G)      --  compute moments of jj correlator G.

        mom2taylor(G)   --  convert moments in G into Taylor coefficients
                            for q2-expansion.

        taylor2mom(G)   --  convert Taylor expansion coefficients G into moments.

        vacpol(G)       --  create a Pade approximant for the subtracted
                            vacuum polarization (PI-hat) from the jj correlator
                            whose moments (or Taylor coefs) are in G.

        a_mu(pihat, Q)  --  compute the contribution to the muon's g-2
                            anomaly from function pihat (usually built by vacpol).

        pade_gvar(f,m,n)--  general-purpose code for determining Pade approximants
                            to a power series whose coefficients are GVars.

        pade_svd(f,m,n) --  general-purpose code for determining Pade approximants
                            for a power series whose coefficients are floats.
                            Uses SVD regularization to stabilize results when
                            the input data are noisy.
    """
    ,
    classifiers = [                     #
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Physics'
        ]

)
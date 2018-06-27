from distutils.command.build_py import build_py as _build_py
from distutils.core import setup

G2TOOLS_VERSION = '1.3.4'

class build_py(_build_py):
    def run(self):
        """ Append version number to g2tools.py """
        with open('src/g2tools.py', 'a') as gfile:
            gfile.write("\n__version__ = '%s'\n" % G2TOOLS_VERSION)
        _build_py.run(self)

setup(name='g2tools',
    version=G2TOOLS_VERSION,
    description='Utilities for muon g-2 analyses in lattice QCD.',
    author='G. Peter Lepage, Cornell University',
    author_email='g.p.lepage@cornell.edu',
    license='GPLv3+',
    packages={''},
    package_dir={'':'src'},
    cmdclass={'build_py': build_py},
    requires=['numpy (>=1.7)', 'gvar (>=7.3)', 'scipy', 'lsqfit'],    # for docutils
    install_requires=['gvar>=7.3', 'numpy>=1.7', 'scipy', 'lsqfit'],  # for pip
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
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Physics'
        ]

)
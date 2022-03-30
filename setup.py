from distutils.core import setup

G2TOOLS_VERSION = open('src/g2tools/_version.py', 'r').readlines()[0].split("'")[1]

# pypi
with open('README.rst', 'r') as file:
    long_description = file.read()

setup(name='g2tools',
    version=G2TOOLS_VERSION,
    description='Utilities for muon g-2 analyses in lattice QCD.',
    author='G. Peter Lepage, Cornell University',
    author_email='g.p.lepage@cornell.edu',
    license='GPLv3+',
    packages={'g2tools'},
    package_dir={'g2tools':'src/g2tools'},
    # cmdclass={'build_py': build_py},
    requires=['numpy (>=1.7)', 'gvar (>=7.3)', 'scipy', 'lsqfit'],    # for docutils
    install_requires=['gvar>=7.3', 'numpy>=1.7', 'scipy', 'lsqfit'],  # for pip
    platforms="Any",
    url="https://github.com/gplepage/g2tools.git",
    long_description=long_description,
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
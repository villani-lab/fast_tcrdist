from setuptools import setup, find_packages
from codecs import open
from os import path
from Cython.Build import cythonize
from distutils.core import setup, Extension

here = path.abspath(path.dirname(__file__))

import numpy
np_include = numpy.get_include()
include_dirs = [np_include, "fast_tcrdist/tcrdist/cython"]

extensions = [Extension('fast_tcrdist.tcrdist.cython.seq_dist', ['fast_tcrdist/tcrdist/cython/seq_dist.pyx']),
Extension("fast_tcrdist.tcrdist.cython.cnwalign", ["fast_tcrdist/tcrdist/cython/cnwalign.c"], include_dirs = include_dirs)]

setup(
    name='fast_tcrdist',

    version='0.01',

    description='Optimized TCRDist calculation for TCR repertoire data analysis',
    long_description='',

    url='',

    author='Neal Smith',
    author_email='nsmith19@mgh.hardvard.edu',

    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'License :: OSI Approved :: Apache Software License',

        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
    ],

    install_requires=["pandas", "numpy", "distance", "cython"],

    keywords='',

    packages=find_packages(),

    entry_points={
        'console_scripts': [
            'sample=sample:main',
        ],
    },
    ext_modules = cythonize(extensions), include_dirs = include_dirs
    # ext_modules = [Extension("TCR_homology_10X.tcrdist.nw_align", ['TCR_homology_10X/tcrdist/nw_align.c'])]
)

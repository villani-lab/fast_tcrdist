from setuptools import setup, find_packages
from codecs import open
from os import path
from Cython.Build import cythonize
from distutils.core import setup, Extension
import numpy

here = path.abspath(path.dirname(__file__))

np_include = numpy.get_include()
include_dirs = [np_include, "fast_tcrdist/tcrdist/cython"]

extensions = [Extension('fast_tcrdist.tcrdist.cython.seq_dist', ['fast_tcrdist/tcrdist/cython/seq_dist.pyx']),
Extension("fast_tcrdist.tcrdist.cython.cnwalign", ["fast_tcrdist/tcrdist/cython/cnwalign.pyx"], include_dirs = include_dirs)]

setup(
    name='fast_tcrdist',

    version='0.0.1',

    description='Optimized TCRDist calculation for TCR repertoire data analysis',
    long_description='',

    url='',

    author='Neal Smith',
    author_email='nsmith19@mgh.harvard.edu',

    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English'
    ],

    install_requires=["pandas", "numpy", "distance", "cython", "hdbscan", "logomaker", "anndata", "scanpy"],

    packages=find_packages(),

    ext_modules = cythonize(extensions), include_dirs = include_dirs
)

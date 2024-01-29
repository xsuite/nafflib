# copyright ############################### #
# This file is part of the Xpart Package.   #
# Copyright (c) CERN, 2021.                 #
# ######################################### #

from setuptools import setup, find_packages, Extension
from pathlib import Path

#######################################
# Prepare list of compiled extensions #
#######################################

extensions = []


#########
# Setup #
#########

version_file = Path(__file__).parent / 'nafflib/_version.py'
dd = {}
with open(version_file.absolute(), 'r') as fp:
    exec(fp.read(), dd)
__version__ = dd['__version__']

setup(
    name='nafflib',
    version=__version__,
    description='nafflib algorith for frequency analysis',
    long_description=('nafflib algorith for frequency analysis'),
    #url='https://xsuite.readthedocs.io/',
    packages=find_packages(),
    ext_modules=extensions,
    include_package_data=True,
    install_requires=[
        'numpy>=1.0',
        'pandas',
        'numba'
        ],
    author='P. Belanger et al.',
    license='Apache 2.0',
    download_url="https://pypi.python.org/pypi/nafflib",
    project_urls={
            "Bug Tracker": "https://github.com/xsuite/xsuite/issues",
            # "Documentation": 'https://xsuite.readthedocs.io/',
            "Source Code": "https://github.com/xsuite/nafflib",
        },
    extras_require={
        'tests': ['pytest'],
        },
    )

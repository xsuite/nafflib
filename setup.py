# copyright ############################### #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from setuptools import setup, find_packages, Extension
from pathlib import Path

#######################################
# Prepare list of compiled extensions #
#######################################

extensions = []

# LOAD REAME as PyPI description
with open("README.md","r") as fh:
    readme = fh.read()

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
    long_description=readme,
    long_description_content_type='text/markdown',
    url='https://github.com/xsuite/nafflib',
    packages=find_packages(),
    ext_modules=extensions,
    include_package_data=True,
    install_requires=[
        'numpy>=1.0',
        'pandas',
        'numba'
        ],
    author='P. Belanger, K. Paraschou et al.',
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

from setuptools import setup, Extension, find_packages
import sys
import numpy as np

if sys.version_info[0] < 3:
    extension_name = "NAFFlib/NAFFlib2_c"
else:
    extension_name = "NAFFlib/NAFFlib_c"

with open("README.md","r") as fh:
    long_description = fh.read()

module = Extension(extension_name,
                   ["NAFFlib/source/brent.c",
                    "NAFFlib/source/fft.c",
                    "NAFFlib/source/frequency.c",
                    "NAFFlib/source/pynafflib.c",
                    "NAFFlib/source/signal_processing.c",
                    "NAFFlib/source/windows.c",
                   ],
                   include_dirs=["NAFFlib/include", np.get_include()]
)

setup(
    name="NAFFlib",
    version="1.0.0",
    author="Konstantinos Paraschou",
    author_email="konstantinos.paraschou@cern.ch",
    description="A Python-wrapped C library which implements the NAFF algorithm",
    long_description=long_description,
    url="https://github.com/PyCOMPLETE/NAFFlib",
    packages=find_packages(),
    license = 'LGPLv2.1',
    keywords = 'frequency analysis naff',
    ext_modules=[module]
)

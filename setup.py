from setuptools import setup, Extension, find_packages
import sys
import numpy as np

if sys.version_info[0] < 3:
    extension_name = "nafflib/nafflib2_c"
else:
    extension_name = "nafflib/nafflib_c"

with open("README.md","r") as fh:
    long_description = fh.read()

module = Extension(extension_name,
                   ["nafflib/source/brent.c",
                    "nafflib/source/fft.c",
                    "nafflib/source/frequency.c",
                    "nafflib/source/pynafflib.c",
                    "nafflib/source/signal_processing.c",
                    "nafflib/source/windows.c",
                   ],
                   include_dirs=["nafflib/include", np.get_include()],
                   extra_compile_args=["-std=c99"]
                  )

setup(
    name="nafflib",
    version="1.0.2",
    author="Konstantinos Paraschou",
    author_email="konstantinos.paraschou@cern.ch",
    description="A Python-wrapped C library which implements the NAFF algorithm",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/xsuite/nafflib",
    packages=find_packages(),
    license = 'LGPLv2.1',
    keywords = 'frequency analysis naff',
    ext_modules=[module]
)

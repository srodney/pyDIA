from setuptools import setup
from setuptools.extension import Extension
import os
import numpy


# extension module(s): 
#fname1 = os.path.join("pydia", "c_functions.c")
#fname2 = os.path.join("pydia", "c_functions_dp.c")
#source_files = [fname1, fname2]
#include_dirs = [numpy.get_include()]
#extensions = [Extension("pydia.c_functions_dp", source_files,
#                        include_dirs=include_dirs)]
extensions = [Extension('pydia/c_functions_dp', sources=['pydia/c_functions_dp.c'])]

setup(
    name='pydia',
    install_requires=['astropy>=0.4.0'],
    author='Michael D. Albrow',
    author_email='',
    version='1.1.0',
    packages=['pydia'],
    zip_safe=False,
    #eager_resources=['pydia/c_functions.so',
    #                 'pydia/c_functions_dp.so'],
    ext_modules=extensions,
)

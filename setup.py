from setuptools import setup
from setuptools.extension import Extension
import os
import shutil
import glob

extensions = [Extension('pydia/c_functions_dp',
                        sources=['pydia/c_functions_dp.c'])]

setup(
    name='pydia',
    install_requires=['astropy>=0.4.0'],
    author='Michael D. Albrow',
    author_email='',
    version='1.1.0',
    packages=['pydia'],
    zip_safe=False,
    ext_modules=extensions,
)

# Catch an unusually-named object file created by cpython on a mac
if not os.path.isfile('pydia/c_functions_dp.so'):
    sofiles = glob.glob("c_functions_dp.*so")
    if len(sofiles)>0:
        sofile = sofiles[0]
        os.rename(sofile, 'pydia/c_functions_dp.so')

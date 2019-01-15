from setuptools import setup
from setuptools.extension import Extension

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

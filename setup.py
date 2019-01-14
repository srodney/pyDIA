# -*- coding: utf-8 -*-

import os
import sys
import setuptools as st
from io import open

# Fix so that the setup.py usage is CWD-independent
SETUPDIR = os.path.abspath(os.path.dirname(__file__))
SETUPDIR = os.path.dirname(__file__)
PKGDIR = os.path.join(SETUPDIR, 'Code')

sys.path.append(PKGDIR)

# Ugh.  Close your eyes!  Making a symlink to Code to make it importable
IMPORTDIR = PKGDIR.replace('Code','pyDIA')
if not os.path.exists(IMPORTDIR):
    os.symlink(PKGDIR, IMPORTDIR)

import pyDIA

reqsfname = os.path.join(SETUPDIR, 'requirements.txt')
reqs = open(reqsfname, 'r', encoding='utf-8').read().strip().splitlines()

longdesc = """pyDIA is a modular python package for performing star detection, 
difference imaging and photometry on astronomical images.
"""

st.setup(
    name="pyDIA",
    version=pyDIA.__version__,
    author=u"",
    author_email="",
    description=("Python Difference Imaging Analysis"),
    license="MIT",
    url="https://github.com/MichaelDAlbrow/pyDIA",
    package_dir = {'': PKGDIR},
    packages = st.find_packages(PKGDIR),
    entry_points = {},
    install_requires=reqs,
    extras_require={},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Natural Language :: English",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Utilities",
        "Topic :: Scientific/Engineering :: Image Processing",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
    ],
    long_description=longdesc,
    zip_safe=True,
)

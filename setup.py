#!/usr/bin/env python3
# Pretty much copied (with a few modifications) from
# https://github.com/astrobc1/pychell/blob/master/setup.py

import setuptools
import os

# Get requirements
theLibFolder = os.path.dirname(os.path.realpath(__file__))
requirementPath = theLibFolder + '/requirements.txt'
install_requires = []
if os.path.isfile(requirementPath):
    with open(requirementPath) as f:
        install_requires = f.read().splitlines()

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="periodogram",
    version="0.1.0",
    author="John Berberian, Jr.",
    author_email="jeb.study@gmail.com",
    description="Compute periodograms of formatted data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=install_requires,
    url="https://github.com/LindirQuenya/Periodogram",
    py_modules=['periodogram', '_periodogram'],
    classifiers=[
        "Programming Language :: Python :: 3.5",
        "Operating System :: Unix",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
    ],
    python_requires='>=3.5'
)

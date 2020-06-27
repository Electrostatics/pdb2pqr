#!/usr/bin/python3

"""
The use of continuum solvation methods such as APBS requires accurate and complete structural data as well as force field parameters such as atomic charges and radii.
Unfortunately, the limiting step in continuum electrostatics calculations is often the addition of missing atomic coordinates to molecular structures from the Protein
Data Bank and the assignment of parameters to these structures. To adds this problem, we have developed PDB2PQR. This software automates many of the common tasks of
preparing structures for continuum solvation calculations as well as many other types of biomolecular structure modeling, analysis, and simulation. These tasks include:

* Adding a limited number of missing heavy (non-hydrogen) atoms to biomolecular structures.
* Estimating titration states and protonating biomolecules in a manner consistent with favorable hydrogen bonding.
* Assigning charge and radius parameters from a variety of force fields.
* Generating "PQR" output compatible with several popular computational modeling and analysis packages.

This service is intended to facilitate the setup and execution of electrostatics calculations for both experts and non-experts and thereby broaden the accessibility of
biomolecular solvation and electrostatics analyses to the biomedical community.
"""

import sys
import setuptools

if sys.version_info[:2] < (3,6):
    raise RuntimeError("Python version >= 3.6 is required.")

with open("README.md", "r") as f:
    long_description = f.read()

MAJOR = 3
MINOR = 0
MICRO = 0
VERSION = "%d.%d.%d" % (MAJOR, MINOR, MICRO)

setuptools.setup(
    name="pdb2pqr",
    version=VERSION,
    author="Jens Erik Nielsen, Nathan Baker, and many others.",
    author_email="nathanandrewbaker@gmail.com",
    description=(
        "Automates many of the common tasks of preparing structures for "
        "continuum solvation calculations as well as many other types of "
        "biomolecular structure modeling, analysis, and simulation."),
    long_description=long_description,
    install_requires=["propka >= 3.2", "pandas >= 1.0", "pytest>=5.4.1"],
    url="https://github.com/Electrostatics/apbs-pdb2pqr/tree/master/pdb2pqr",
    packages=setuptools.find_packages(
        exclude=["pdb2pka", "*.pdb2pka", "pdb2pka.*", "*.pdb2pka.*"]),
    package_data={"pdb2pqr": ["dat/*.xml", "dat/*.DAT", "dat/*.names"]},
    python_requires=">=3.5",
    license="BSD",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: Linux",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemestry"
    ],
    keywords="science chemistry molecular biology",
    entry_points={"console_scripts": "pdb2pqr30=pdb2pqr.main:main"}
)

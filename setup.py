from setuptools import setup
import os
from os import path
import sys


this_directory = path.abspath(path.dirname(__file__))
sys.path.append(this_directory)

import mdpertool

with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mdpertool', '_version.py')) as version_file:
    exec(version_file.read())

install_requires = [
    'python==3.11',
    'numpy',
    'pyvis',
    'matplotlib',
    'pandas',
    'networkx',
    'pyopengl',
    'biopython',
    'prody',
    'pyqtwebengine',
    'pyyaml',
    'pystache',  
]

# Conda-forge ile yüklenen paketler
conda_packages = [
    "mdtraj",
    "openmm",
    "pymol-open-source",
    "pyside2",
    "mdanalysis",
    "pyqtgraph",
    "parmed",
]


setup(
    name='mdpertool',
    version=__version__,
    description=__description__,
    long_description=__long_description__,
    url=__url__,
    author=__author__,
    author_email=__author_email__,
    license='MIT',
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "License :: OSI Approved :: Academic Free License (AFL)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        'Programming Language :: Python',
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Software Development :: User Interfaces",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    packages=['mdpertool'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'mdpertool=mdpertool.ui_main:run_mdpertool',
        ],
    },

    extras_require={
        'conda': conda_packages
    }
)
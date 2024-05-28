from setuptools import setup, find_packages
import os
import os.path as op
import sys


with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mdpertool', '_version.py')) as version_file:
    exec(version_file.read())


install_requires = [
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
    "cudatoolkit",
    "mdanalysis",
    "pyqtgraph",
    "parmed",
]

tests_require = [
    'pytest',
    'pytest-cov',
]

needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
setup_requires = ['pytest-runner'] if needs_pytest else []

def find_package_data():
    files = []
    for root, dirnames, filenames in os.walk('mdpertool'):
        for fn in filenames:
            files.append(os.path.relpath(os.path.join(root, fn), 'mdpertool'))

    return files

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
    packages=find_packages(),
    package_data={'mdpertool': find_package_data()},
    platform='any',
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'mdpertool=mdpertool.ui_main:run_mdpertool',
        ],
    },

    install_requires=install_requires,
    tests_require=tests_require,
    setup_requires=setup_requires,
    extras_require={
        'conda': conda_packages
    }
)

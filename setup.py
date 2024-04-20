from setuptools import setup, find_packages
import os
import sys

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mdpertool', '_version.py')) as version_file:
    exec(version_file.read())

install_requires = [
    'numpy',
    'biopython',
    'matplotlib',
    'mdtraj',
    'networkx',
    'pandas',
    'parmed',
    'prody',
    'pyopengl',
    'pyqtgraph',
    'pyqtwebengine',
    'pystache',
    'pyvis',
    ]

# Conda-forge ile yüklenen paketleri, dependency_links kullanarak belirtebilirsiniz.
dependency_links=[
    "git+https://github.com/conda-forge/mdtraj-feedstock",
    "git+https://github.com/conda-forge/openmm-feedstock",
    "git+https://github.com/conda-forge/pymol-open-source-feedstock",
    "git+https://github.com/conda-forge/pyside2-feedstock",
    "git+https://github.com/conda-forge/cudatoolkit-feedstock",
    "git+https://github.com/conda-forge/mdanalysis-feedstock"
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
        "Development Status :: 2 - Pre-Alpha",
        "Development Status :: 4 - Beta",
        "Environment :: GPU",
        "Environment :: GPU :: NVIDIA CUDA",
        "Environment :: MacOS X",
        "License :: OSI Approved :: MIT License",
        "License :: OSI Approved :: Academic Free License (AFL)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: Unix",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Software Development :: Libraries",
        "Topic :: Software Development :: User Interfaces",
        "Topic :: Software Development :: Widget Sets",
    ],
    packages=find_packages(),
    package_data={'mdpertool': find_package_data()},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'mdpertool=mdpertool.ui_main:run_mdpertool',
        ],
    },

    install_requires=install_requires,
    tests_require=tests_require,
    setup_requires=setup_requires,
    dependency_links=dependency_links
)

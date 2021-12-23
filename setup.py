import os

from setuptools import setup, find_packages
import os.path as op
import sys

with open(op.join(op.dirname(op.realpath(__file__)), 'mdpertool', '_version.py')) as version_file:
    exec(version_file.read())
#file_path = op.join(op.dirname(op.realpath(__file__)), 'mdpertool', '_version.py')
#os.system(f'py {file_path}')

with open(op.join(op.dirname(op.realpath(__file__)), 'README.md')) as readme_file:
    readme = readme_file.read()

with open(op.join(op.dirname(op.realpath(__file__)), 'CHANGELOG.md')) as changelog_file:
    changelog = changelog_file.read()

with open(op.join(op.dirname(op.realpath(__file__)), 'CITATION.md')) as citation_file:
    citation = citation_file.read()

desc = readme + '\n\n' + changelog + '\n\n' + citation
try:
    import pypandoc

    long_description = pypandoc.convert_text(desc, 'rst', format='md')
    with open(op.join(op.dirname(op.realpath(__file__)), 'README.rst'), 'w') as rst_readme:
        rst_readme.write(long_description)
except (ImportError, OSError, IOError):
    long_description = desc

install_requires = [
    'openmm >=7.6',
    'numpy >=1.20',
    'pyvis >=0.19',
    'pymol-open-source >=2.5',
    'pyyaml',
    'matplotlib >=3.4',
    'pyopengl >=3.1',
    'pdbfixer >=1.8',
    'mdtraj >=1.9.5',
    'networkx >=2.6',
    'parmed >=3.2',
    'prody >=2.0',
    'pyside2 >=5.13',
    'pystache >=0.5',
    'pyqtgraph >=0.12'
]

tests_require = [
    'pytest>=3.0.0',
    'pytest-cov>=2.3.1',
]

needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
setup_requires = ['pytest-runner'] if needs_pytest else []

# github_token = os.environ['GITHUB_TOKEN']

setup(
    name='mdpertool',
    version=__version__,
    description='A Software Tool for Investigation of Allosteric Communication within Protein Structures via Energy '
                'Dissipation in Molecular Dynamics Simulations',
    long_description=long_description,
    url='https://github.com/bio-otto/MDPerTool_GUI',
    author=__author__,
    author_email='halil.ibrahim.oozdemir@gmail.com',
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
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Software Development :: Libraries",
        "Topic :: Software Development :: User Interfaces",
        "Topic :: Software Development :: Widget Sets",
    ],
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'processrcmfolder=mdpertool.ui_main:ui_main.py',
        ],
    },
    install_requires=install_requires,
    tests_require=tests_require,
    setup_requires=setup_requires,
    # dependency_links=[
    #     "openmm @ git+https://github.com/openmm/openmm",
    #     #               "parmed @ git+https://github.com/ParmEd/ParmEd",
    #     # 'pdbfixer @ git+https://github.com/openmm/pdbfixer@master',
    #     # "pdbfixer @ git+https://github.com/openmm/pdbfixer",
    #     # "pymol-open-source @ git+https://github.com/schrodinger/pymol-open-source"
    # ]
    # dependency_links=['https://github.com/openmm/pdbfixer/master#egg=pdbfixer-1.8.1',
    #                   'https://github.com/schrodinger/pymol-open-source/master#egg=pymol-2.5.0',
    #                   'https://github.com/openmm/openmm/master#egg=openmm-7.6.0',
    #                   'https://github.com/ParmEd/ParmEd/master#egg=parmed-3.4.3',
    #                   'https://github.com/mdtraj/mdtraj/master#egg=mdtraj-1.9.7']
)

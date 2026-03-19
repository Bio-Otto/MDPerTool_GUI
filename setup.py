from pathlib import Path

from setuptools import find_namespace_packages, setup


this_directory = Path(__file__).resolve().parent
version_ns = {}
exec((this_directory / 'mdpertool' / '_version.py').read_text(encoding='utf-8'), version_ns)
long_description = (this_directory / 'README.md').read_text(encoding='utf-8')


# Conda-forge ile yüklenen paketler

conda_packages = [

    "mdtraj<=1.9.9",
    "openmm>=8.0,<8.2",
    "pymol-open-source>=2.5",
    "pyside2>=5.15,<6",
    "mdanalysis",
    "pyqtgraph",
    "parmed",
]


setup(

    name='mdpertool',
    version=version_ns['__version__'],
    description=version_ns['__description__'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    url=version_ns['__url__'],
    author=version_ns['__author__'],
    author_email=version_ns['__author_email__'],
    license='MIT',
    python_requires='>=3.9,<3.10',
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

    packages=find_namespace_packages(include=['mdpertool', 'mdpertool.*']),


    install_requires=[
        'prody',
    ],


    include_package_data=True,

    package_data={

        'mdpertool': [
            'Download/*',
            'analysis/*',
            'fonts/*',
            'gui/*',
            'no_gui/*',
            'src/*'

        ]

    },


    entry_points={

        'console_scripts': [

            'mdpertool=mdpertool.launcher:run_mdpertool',

        ],

    },


    extras_require={
        'conda': conda_packages

    }

)
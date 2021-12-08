# try:
#     from setuptools import setup
# except ImportError:
#     from distutils.core import setup
#
# setup(name='mdpertool',
#       version='0.1',
#       description='Perturbation based Allosteric Pathway Finder',
#       url="https://github.com/Bio-Otto/MDPERTOOL_v01/tree/gui_development",
#       author='Halil Ibrahim Ozdemir',
#       author_email='halilibrahim307@marun.edu.tr',
#       license='MIT',
#       install_requires=[
#           'pyside2',
#           'mdanalysis',
#           'mdtraj',
#           'matplotlib',
#           'pandas',
#           'numpy',
#           'pyqtgraph',
#       ],
#       zip_safe=False)

from setuptools import setup, find_packages

description = "Perturbation based Allosteric Pathway Finder"

long_description = """A Software Tool for Investigation of Allosteric Communication within Protein Structures via Energy Dissipation in Molecular Dynamics Simulations"""

author = "Halil Ibrahim Ozdemir"

maintainer = author

maintainer_email = "halilibrahim307@marun.edu.tr"

url = "https://github.com/Bio-Otto/MDPERTOOL_v01/tree/gui_development"

download_url = url

platforms = ["Linux", "Windows", "Mac"]

keywords = ["Perturbation", "Energy", "Dissipation", "Allostery", "OpenMM", "Molecular Dynamics"]

install_requires = ["pyqtgraph", "parmed", "mdtraj", "matplotlib", "pandas", "numpy", "openmm", "pystache",
                    "pyside2==5.13.2", "pymol-open-source", "networkx", "pdbfixer", "prody", "pyvis", "pyopengl"]



classifiers = [
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
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Visualization",
    "Topic :: Software Development :: Libraries",
    "Topic :: Software Development :: User Interfaces",
    "Topic :: Software Development :: Widget Sets",
]


setup(
    name="MDPERTOOL_v01",
    version="0.1",
    description=description,
    long_description=long_description,
    author=author,
    maintainer=maintainer,
    maintainer_email=maintainer_email,
    url=url,
    download_url=download_url,
    platforms=platforms,
    license="MIT",
    keywords=keywords,
    packages=find_packages(exclude=["tests"]),
    classifiers=classifiers,
    include_package_data=True,
    test_suite="tests",
    python_requires=">=3.7",
    install_requires=install_requires,
    zip_safe=False
)

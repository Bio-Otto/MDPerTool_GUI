# Installation

MDPerTool can be easily installed using `conda`. We highly recommend using a fresh virtual environment.

## Requirements
- Anaconda or Miniconda
- Windows, macOS, or Linux

## Installation via Conda

```bash
conda create -n mdpertool_env python=3.9
conda activate mdpertool_env
conda install -c conda-forge -c defaults openmm mdtraj
conda install -c bio-otto mdpertool
```

Verify the installation:
```bash
mdpertool --help
```
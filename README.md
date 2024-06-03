# MDPERTOOL: Perturbation-based Allosteric Pathway Finder

<a href="http://compbio.bioe.eng.marmara.edu.tr/" target="_parent">
<img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/master/mdpertool/gui/icons/big_icons/logo_for_contacts.png" style="position:fixed; bottom:-15px; width: 150px; height: 120px;" align="right"/>
</a>



<div align="center">
	<a href="https://anaconda.org/bio-otto/mdpertool"><img src="https://img.shields.io/conda/v/bio-otto/mdpertool.svg" alt="Conda Build"/></a>
    <img src="https://anaconda.org/bio-otto/mdpertool/badges/platforms.svg" alt="Platforms"/>
    <img src="https://anaconda.org/bio-otto/mdpertool/badges/downloads.svg" alt="Download Counter"/>
    <img src="https://anaconda.org/bio-otto/mdpertool/badges/license.svg" alt="License"/>
    <a href="https://github.com/Bio-Otto/MDPerTool_GUI/contributors"><img alt="GitHub contributors" src="https://img.shields.io/github/contributors/Bio-Otto/MDPerTool_GUI?color=2b9348"></a>
    <a href="http://compbio.bioe.eng.marmara.edu.tr/"><img src="https://img.shields.io/static/v1?label=&labelColor=505050&message=ozbek-lab&color=%230076D6&style=flat&logo=google-chrome&logoColor=%230076D6" alt="website"/></a>
    <a href="https://twitter.com/LabOzbek"><img src="https://img.shields.io/twitter/follow/pemoshh.svg?style=social" alt="Twitter Follow"/></a>
    <a href="https://github.com/Bio-Otto/MDPerTool_GUI/stargazers"><img src="https://img.shields.io/github/stars/Bio-Otto/Bio-Otto" alt="Stars Badge"/></a>
    <a href="https://github.com/Bio-Otto/MDPerTool_GUI/members"><img src="https://img.shields.io/github/forks/Bio-Otto/Bio-Otto" alt="Forks Badge"/></a>
   
</div>

<i>A Software Tool for Investigation of Allosteric Communication within Protein Structures via Energy Dissipation in Molecular Dynamics Simulations</i>

![MDPerTool Animation](https://github.com/Bio-Otto/MDPerTool_GUI/blob/master/mdpertool/gui/icons/MDPerTool.gif)

---

## 📥 Installation

To install MDPERTOOL using Conda, simply run the following commands:

```sh
conda install bio-otto::mdpertool
```
This will install MDPERTOOL and all of its dependencies. Once the installation is completed, you can access the Command Line Interface by typing 'mdpertool' in your terminal.

---

### Manual Installation

To manually install MDPERTOOL, follow the instructions below:

#### For Windows Users
```cmd
git clone https://github.com/Bio-Otto/MDPerTool_GUI.git
cd MDPerTool_GUI
conda env create -f envs/winx64_env.yml
conda activate mdpertool
cd mdpertool
python ui_main.py
```
<p align="right">
	<img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/master/mdpertool/gui/icons/windows-logo.png" width="60" title="Available on Windows">
</p>


#### For Linux Users
```sh
git clone https://github.com/Bio-Otto/MDPerTool_GUI.git
cd MDPerTool_GUI
conda env create -f envs/linux64_env.yml
conda activate mdpertool
cd mdpertool
python ui_main.py
```
<p align="right">
    <img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/master/mdpertool/gui/icons/orange-logo-linux.png" width="35" title="Available on Ubuntu 20.10">
</p>

---

## 🖥️ Terminal Usage

For terminal usage, run:

```sh
python no_gui.py -p pdb_file -wdcd True -pert_res 'SER345' -speed_factor 4
```


## 📦 Dependencies

* __OpenMM__ - A high performance toolkit for molecular simulation. 
* __Networkx__ - NetworkX is a Python package for the creation, manipulation, and study of the structure, dynamics, and functions of complex networks.

* __Pandas__ - pandas is a fast, powerful, flexible and easy to use open source data analysis and manipulation tool
* __Numpy__ -  The fundamental package for scientific computing with Python 
* __Pymol (Open Source)__ - PyMOL is a user-sponsored molecular visualization system on an open-source foundation, maintained and distributed by Schrödinger.
* __Matplotlib__ - Matplotlib is a comprehensive library for creating static, animated, and interactive visualizations in Python.
* __Pyqtgraph__ - Scientific Graphics and GUI Library for Python
* __PySide2__ - Python bindings for the Qt cross-platform application and UI framework
* __ProDy__ - Protein Dynamics and Sequence Analysis
* __Parmed__ - Parameter/topology editor and molecular simulator

And of course MDPERTOOL itself is an open source public repository.

#### Also you can check full functional parameters with: 

```sh
python no_gui.py -h
```

For Molecular Dynamic Simulation production

```sh
python no_gui.py -p <pdb file> -pff <charmm36> -wff <tip5p> -wdcd <True>
```

MDPerTool applies Energy Dissipation Concept using OpenMM Molecular Dynamic Toolkit, which also supports the CUDA platform.


## ✨ Features


  - Automated topology builder
  - Residue base energy decomposition
  - DCD and XTC file format support
  - Easy to installation using Conda 
  

### Features of MDPERTOOL v0.0.1

This version of mdpertool includes the following features. 

<table align="center">
    <thead>
        <tr>
            <th align="left">Feature</th>
            <th align="center">Documentation</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td align="left">Each Residue Decomposition</td>
            <td align="center">[plugins/dropbox/README.md][MDPERTOOL]</td>
        </tr>
        <tr>
            <td align="left">Molecular Dynamic Simulation</td>
            <td align="center">[plugins/github/README.md][MDPERTOOL]</td>
        </tr>
        <tr>
            <td align="left">Energy Dissipation Network</td>
            <td align="center">[plugins/googledrive/README.md][MDPERTOOL]</td>
        </tr>
        <tr>
            <td align="left">Free Energy Calculations</td>
            <td align="center">[plugins/onedrive/README.md][MDPERTOOL]</td>
        </tr>
    </tbody>

</table><p></p></div>

---

## 🚀 Development

Want to contribute? Get a branch and join us to make MDPERTOOL even greater!


## 📋 Todos

 - Write MORE Tests
 - Add Light Mode


## 📝 License

MDPERTOOL is licensed under the MIT License. See the LICENSE file for more details.

Free Software, Hell Yeah!
<p align="center">
    <a href="http://compbio.bioe.eng.marmara.edu.tr/" target="_parent">
    <img src="https://github.com/Bio-Otto/Example_MD_Scripts/blob/master/PoweredByOzbekLab.png" width="250" height="40" /></a>
    <img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/master/mdpertool/gui/icons/orange-logo-linux.png" width="40" title="Available on Ubuntu 20.10">
    <img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/master/mdpertool/gui/icons/windows-logo.png" width="40" title="Available on Windows">
</p>
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
This will install MDPERTOOL and all of its dependencies. Once the installation is completed, you can access the Command Line Interface by typing mdpertool cli in your terminal or the Graphical User Interface by typing mdpertool gui.

---

### Manual Installation

To manually install MDPERTOOL, follow the instructions below:

#### For Windows Users
```sh
git clone https://github.com/Bio-Otto/MDPerTool_GUI.git
cd MDPerTool_GUI
conda env create -f envs/winx64_env.yml
conda activate mdpertool
cd mdpertool
python ui_main.py
```
<p align="right">
	<img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/gui_development/mdpertool/gui/icons/windows-logo.png" width="70" title="Available on Windows">
</p>


#### For Linux and macOS Users
```sh
git clone https://github.com/Bio-Otto/MDPerTool_GUI.git
cd MDPerTool_GUI
conda env create -f envs/linux64_env.yml
conda activate mdpertool
cd mdpertool
python ui_main.py
```
<p align="right">
    <img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/gui_development/mdpertool/gui/icons/orange-logo-linux.png" width="35" title="Available on Ubuntu">
    <img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/gui_development/mdpertool/gui/icons/macOS-logo.png" width="50" title="Available on macOS">
</p>

---

## 🖥️ Usage

### Graphical User Interface (GUI)

To launch the GUI, simply run:

```sh
mdpertool gui
```

### Command Line Interface (CLI)

To run the CLI mode, you need to provide additional arguments. The basic command structure is:

```sh
mdpertool cli -p <topology_file> [optional arguments]
```

#### Example Usage

Running in CLI Mode
```sh
mdpertool cli -p path/to/topology.pdb -pff amber99sb -wff spce
```
This command will run MDPERTOOL in command line mode with the specified topology file, protein forcefield, and water forcefield.

#### Full CLI Options

For a complete list of CLI options, you can refer to the help command:
```sh
mdpertool cli -h
```
This will display all the available options and their descriptions for running MDPERTOOL in CLI mode.


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
    <img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/gui_development/mdpertool/gui/icons/orange-logo-linux.png" width="40" title="Available on Ubuntu">
    <img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/gui_development/mdpertool/gui/icons/windows-logo.png" width="70" title="Available on Windows">
    <img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/gui_development/mdpertool/gui/icons/macOS-logo.png" width="50" title="Available on macOS">
</p>
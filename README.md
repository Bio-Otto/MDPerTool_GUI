
<h1 align="left"><span style="" text-align="left">MDPERTOOL: Perturbation based Allosteric Pathway Finder</span>
<a href="http://compbio.bioe.eng.marmara.edu.tr/" target="_parent">
<img src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/gui_development/icons/big_icons/logo_for_contacts.png" style="position:fixed; bottom:0px; width: 150px; height: 120px;" width="150" height="120" align="right"/>
</a>

</h1>
<div align="center">
<img src="https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg" alt="Awesome Badge"/>
<a href="http://compbio.bioe.eng.marmara.edu.tr/"><img src="https://img.shields.io/static/v1?label=&labelColor=505050&message=ozbek-lab&color=%230076D6&style=flat&logo=google-chrome&logoColor=%230076D6" alt="website"/></a>
<!-- <img src="http://hits.dwyl.com/abhisheknaiidu/awesome-github-profile-readme.svg" alt="Hits Badge"/> -->
<img src="https://img.shields.io/static/v1?label=%F0%9F%8C%9F&message=If%20Useful&style=style=flat&color=BC4E99" alt="Star Badge"/>
<a href="https://twitter.com/LabOzbek" ><img src="https://img.shields.io/twitter/follow/pemoshh.svg?style=social" /> </a>
<br>

<i>A Software Tool for Investigation of Allosteric Communication within Protein Structures via Energy Dissipation in Molecular Dynamics Simulations</i>

<a href="https://github.com/Bio-Otto/MDPerTool_GUI/stargazers"><img src="https://img.shields.io/github/stars/Bio-Otto/Bio-Otto" alt="Stars Badge"/></a>
<a href="https://github.com/Bio-Otto/MDPerTool_GUI/members"><img src="https://img.shields.io/github/forks/Bio-Otto/Bio-Otto" alt="Forks Badge"/></a>
<a href="https://github.com/Bio-Otto/MDPerTool_GUI/pulls"><img src="https://img.shields.io/github/issues-pr/Bio-Otto/Bio-Otto" alt="Pull Requests Badge"/></a>
<a href="https://github.com/Bio-Otto/MDPerTool_GUI/issues"><img src="https://img.shields.io/github/issues/Bio-Otto/Bio-Otto" alt="Issues Badge"/></a>
<a href="https://github.com/Bio-Otto/MDPerTool_GUI/contributors"><img alt="GitHub contributors" src="https://img.shields.io/github/contributors/Bio-Otto/MDPerTool_GUI?color=2b9348"></a>
<a href="https://github.com/Bio-Otto/MDPerTool_GUI/blob/master/LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License Badge"/></a>
<a href="https://github.com/Bio-Otto/MDPerTool_GUI/blob/master/LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License Badge"/></a>
![MDPERTOOL v0.1 ff69b4](https://img.shields.io/badge/MDPERTOOL-v0.1-<ff69b4>)

<img alt="Awesome GitHub Profile Readme" src="https://github.com/Bio-Otto/MDPerTool_GUI/blob/gui_development/icons/MDPerTool.gif"> </img>



<p align="center">
    <a href="http://compbio.bioe.eng.marmara.edu.tr/" target="_parent">
    <img src="https://github.com/Bio-Otto/Example_MD_Scripts/blob/master/PoweredByOzbekLab.png" style="position:fixed; bottom:0px; width: 250px; height: 40px;" width="250" height="40" /></a>
    <img src="https://github.com/Bio-Otto/MDPERTOOL_v01/blob/ubuntu_gui_development/icons/orange-logo-linux.png" width="40" title="Available on Ubuntu 20.10">
    <img src="https://github.com/Bio-Otto/MDPERTOOL_v01/blob/ubuntu_gui_development/icons/windows-logo.png" width="40" title="Available on Windows">
</p>

### Installation

Open your favorite Terminal and run these commands.


```diff
$ git clone https://github.com/Bio-Otto/MDPerTool_GUI.git
$ cd cd MDPERTOOL_v01
$ conda env create -f envs/env.yml 
$ conda activate mdpertool
$ python ui_main.py
```

### For terminal usage just use; 

```sh
$ python no_gui.py -p pdb_file -wdcd True -pert_res 'SER345' -speed_factor 4
```

* __-p__  -->  Give the absolute path of your pdb file. 
* __-wdcd__  -->  The program defaultly will use dcd reporting. But you can exchange it with XTC file format using ```-wdcd False -wxtc True``` argument
* __-pert_res__  -->  You must list the residue or residues you want to perturbed.
* __-speed_factor__  -->  Indicate how many times you want to increase the velocity of the residue atoms you want to perturbed.


## Then run it.

### Dependencies

MDPERTOOL uses a number of open source projects to work properly:

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

And of course MDPERTOOL v0.1 itself is open source with a [public repository][MDPERTOOL] on GitHub.

### Also you can check full functional parameters with typing 

```sh
$ python no_gui.py -h
```

For production Molecular Dynamic Simulation just type...

```sh
$ no_gui.py -p <pdb file> -pff <charmm36> -wff <tip5p> -wdcd <True>
```

The Program applying Energy Dissipation Concept using powerfull OpenMM Molecular Dynamic Toolkit, which also supports the Cuda platform. Each residual energy calculation required for the concept can be calculated using OpenMMs flexible and useful infrastructure.In addition, you can use the package only for energy decomposition. For this, it will be sufficient to specify a XTC or a DCD file in the script.

# New Features!

  - Automated Topology Builder
  - Residue Decomposition
  - DCD and XTC file format support 
  

### Features of MDPERTOOL v0.1

MDPERTOOL is currently extended with the following features. Instructions on how to use them in your own works are linked below.

<table align="center">
    <thead>
        <tr>
            <th align="left">Feature</th>
            <th align="center">README</th>
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


### Development

Want to contribute?
Get branch and Join us to make MDPERTOOL great!

### Todos

 - Write MORE Tests
 - Add Night Mode

License
----

MIT


**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


[MDPERTOOL]: <https://github.com/Bio-Otto/MDPERTOOL_v01>

# Energy Dissipation Concept with MDPERTOOL v0.1

[![Powered by |Ozbek' Lab](https://github.com/Bio-Otto/Example_MD_Scripts/blob/master/PoweredByOzbekLab.png)](http://compbio.bioe.eng.marmara.edu.tr/)

A Software Tool for Investigation of Allosteric Communication within Protein Structures via Energy Dissipation in Molecular Dynamics Simulations

For terminal usage just use **_python no_gui.py -p pdb_file -wdcd True -pert_res 'SER345' -speed_factor 4_**

* **_-p_**    --> Give the absolute path of your pdb file. 
* **_-wdcd_**    --> The program defaultly will use dcd reporting. But you can exchange it with XTC file format using **_-wdcd False -wxtc True_** argument
* **_-pert_res_**   --> You must list the residue or residues you want to perturbed.
* **_-speed_factor_**   --> Indicate how many times you want to increase the velocity of the residue atoms you want to perturbed.


## Then run it.

### Dependicies

MDPERTOOL uses a number of open source projects to work properly:

* __OpenMM__ - A high performance toolkit for molecular simulation. 
* __Networkx__ - NetworkX is a Python package for the creation, manipulation, and study of the structure, dynamics, and functions of complex networks.

* __Pandas__ - pandas is a fast, powerful, flexible and easy to use open source data analysis and manipulation tool
* __Numpy__ -  The fundamental package for scientific computing with Python 
* __Pymol (Open Source)__ - PyMOL is a user-sponsored molecular visualization system on an open-source foundation, maintained and distributed by Schrödinger.
* __Matplotlib__ - Matplotlib is a comprehensive library for creating static, animated, and interactive visualizations in Python.
* __Mayavi__ - 3D scientific data visualization and plotting in Python
* __Pyqt5__ - Python bindings for the Qt cross platform application toolkit

And of course MDPERTOOL v0.1 itself is open source with a [public repository][MDPERTOOL] on GitHub.

### Todos

 - Write MORE Tests
 - Add Night Mode

License
----

MIT


**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


[MDPERTOOL]: <https://github.com/Bio-Otto/MDPERTOOL_v01>

### Also you can check full functional parameters with typing <**_python no_gui.py -h_**>


The Program applying Energy Dissipation Concept using powerfull OpenMM Molecular Dynamic Toolkit, which also supports the Cuda platform. Each residual energy calculation required for the concept can be calculated using OpenMMs flexible and useful infrastructure.In addition, you can use the package only for energy decomposition. For this, it will be sufficient to specify a XTC or a DCD file in the script.



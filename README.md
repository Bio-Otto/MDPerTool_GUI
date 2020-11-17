# Energy Dissipation Concept with MDPERTOOL v0.1

A Software Tool for Investigation of Allosteric Communication within Protein Structures via Energy Dissipation in Molecular Dynamics Simulations

For terminal usage just use **_python no_gui.py -p pdb_file -wdcd True -pert_res 'SER345' -speed_factor 4_**

* **_pdb_file_** --> Give the absolute path of your pdb file. 
* **_-wdcd:_** --> The program defaultly will use dcd reporting. But you can exchange it with XTC file format using **_-wdcd False -wxtc True_** argument
* **_-pert_res:_** --> You must list the residue or residues you want to perturbed.
* **_-speed_factor:_** --> Indicate how many times you want to increase the velocity of the residue atoms you want to perturbed.


## Then run it.

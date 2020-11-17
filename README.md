# Energy Dissipation Concept with MDPERTOOL v0.1

A Software Tool for Investigation of Allosteric Communication within Protein Structures via Energy Dissipation in Molecular Dynamics Simulations

For terminal usage just use **'<python no_gui.py -p pdb_file -wdcd True -pert_res 'SER345' -speed_factor 4>'**

*-wdcd:* --> The program defaultly will use dcd reporting. But you can exchange it with XTC file format using -wdcd False -wxtc True''' argument
*-pert_res:* --> You must list the residue or residues you want to perturbed.
*-speed_factor:* --> Indicate how many times you want to increase the velocity of the residue atoms you want to perturbed.


                        

Give the absolute path of your pdb file (indicated as *pdb_file*). Then run it.

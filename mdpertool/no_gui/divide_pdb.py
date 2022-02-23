import MDAnalysis as mda

topology = "gast_10.pdb"
multi_frame_pdb = "gas_output.pdb"

u = mda.Universe(topology, multi_frame_pdb)


protein = u.select_atoms("protein")
with mda.Writer("protein.pdb", protein.n_atoms) as W:
    for k, ts in enumerate(u.trajectory):
        if k == 10:
            W.write(protein)


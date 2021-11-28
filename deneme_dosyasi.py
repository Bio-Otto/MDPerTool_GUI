import mdtraj as md
s = 'C:\\Users\\HIbrahim\\Desktop\\out\\without_energy_perturbation_trajectory.dcd'

traj1 = md.join(md.iterload(s, chunk=100, stride=1, atom_indices=None, top='C:\\Users\\HIbrahim\\Desktop\\out\\last.pdb'))

print(traj1)
import mdtraj as md
import numpy as np
#Load and image to remove translational and rotational effects.
traj = md.load('seg.nc', top='mod.prmtop')
ref = md.load('diala.pdb')
traj = traj.image_molecules()
bb = [1, 4, 5, 6, 8, 10, 14, 15, 16, 19]
traj.superpose(ref, atom_indices=bb)
new_traj = traj.atom_slice(bb)
new_traj_ff = new_traj[-1]
coords = np.array(new_traj_ff.xyz)
np.save('coordinates.npy', coords)

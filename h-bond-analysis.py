from Bio.PDB import *
import numpy as np

# Define the PDB file and DSSP file for the MD trajectory
pdb_file = "trajectory.pdb"
dssp_file = "trajectory.dssp"

# Define the chain IDs for the protein and ligand
protein_chain = "A"
ligand_chain = "B"

# Define the atom names for the hydrogen bond donor and acceptor
donor_atom = "OD1"
acceptor_atom = "ND2"

# Define the cutoff distances and angles for the hydrogen bond
distance_cutoff = 3.5
angle_cutoff = 120

# Load the PDB file and DSSP file
parser = PDBParser()
structure = parser.get_structure("trajectory", pdb_file)
model = structure[0]
dssp = DSSP(model, dssp_file)

# Select the protein and ligand chains
protein = model[protein_chain]
ligand = model[ligand_chain]

# Define a list to store the hydrogen bond distances and angles
hbond_distances = []
hbond_angles = []

# Iterate over the MD trajectory
for residue in protein:
    # Only consider residues that are within 5 Angstroms of the ligand
    if residue.get_resname() in ["ARG", "ASN", "ASP", "GLN", "GLU", "HIS", "LYS", "SER", "THR", "TYR", "TRP"]:
        for atom in residue:
            # Only consider hydrogen bond donor atoms
            if atom.get_name() == donor_atom:
                donor_coord = atom.get_coord()
                for altloc in atom.child_dict:
                    for hbond in altloc.child_dict[acceptor_atom].get_bonds():
                        # Only consider hydrogen bond acceptor atoms
                        if hbond.child.get_name() == acceptor_atom:
                            acceptor_coord = hbond.child.get_coord()
                            distance = np.linalg.norm(donor_coord - acceptor_coord)
                            if distance <= distance_cutoff:
                                # Calculate the angle between the hydrogen bond donor, acceptor, and their neighboring atoms
                                neighbor_coords = []
                                for neighbor in hbond.child.get_parent().child_list:
                                    if neighbor.get_name() != acceptor_atom:
                                        neighbor_coords.append(neighbor.get_coord())
                                if len(neighbor_coords) == 2:
                                    angle = calc_angle(neighbor_coords[0], acceptor_coord, neighbor_coords[1])
                                    if angle >= angle_cutoff:
                                        hbond_distances.append(distance)
                                        hbond_angles.append(angle)

# Print the number of hydrogen bonds and their average distance and angle
num_hbonds = len(hbond_distances)
if num_hbonds > 0:
    avg_distance = sum(hbond_distances) / num_hbonds
    avg_angle = sum(hbond_angles) / num_hbonds
    print("Found", num_hbonds, "hydrogen bonds with an average distance of", avg_distance, "Angstroms and an average angle of", avg_angle, "degrees.")
else:
    print("No hydrogen bonds found.")

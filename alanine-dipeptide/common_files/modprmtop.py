import numpy as np
from decimal import Decimal
import sys

# Read input arguments
lamb = float(sys.argv[1])
lookup_str = sys.argv[2]
input_file = sys.argv[3]
output_file = sys.argv[4]
# Number of columns in the section
num_columns = 5

# Read all the lines from the input file
with open(input_file) as file:
    lines = file.readlines()

sent = False
ind_counter = 0
Lennard_jones_count = 0

# Define indices that should be scaled

# Read LJ index pairs from the topology file

# Loop through all the lines in the file to find the section
sent = False
ind_counter = 0
scale_entries = set()

if lookup_str == 'CHARGE':
    try:
        atom_count = int(sys.argv[5])
    except:
        print("Need the number of atoms in the system that are not water!")
        exit()
if lookup_str in ['LENNARD_JONES_ACOEF', 'LENNARD_JONES_BCOEF']:
    try:    
        last_solute_LJTypes = int(sys.argv[5])
        last_solvent_LJTypes = int(sys.argv[6])
        scale_indices = set(range(1, last_solute_LJTypes+1))  #  should always be scaled
        exclude_indices = set(range(last_solute_LJTypes+1, last_solvent_LJTypes+1))# should NOT be scaled exclusively
        lj_index_pairs = []

        in_nb_index_section = False
        for line in lines:
            if "NONBONDED_PARM_INDEX" in line:
                in_nb_index_section = True
                continue
            if in_nb_index_section and "FLAG" in line:
                break
            if in_nb_index_section:
                clean_values = [val for val in line.split() if val.replace("-", "").isdigit()]  # Filter valid numbers
                lj_index_pairs.extend(map(int, clean_values))

        # Identify entries that should be scaled
        for idx, pair_val in enumerate(lj_index_pairs):
            i = (idx // last_solvent_LJTypes) + 1  # First atom type
            j = (idx % last_solvent_LJTypes) + 1  # Second atom type

            if (i in scale_indices or j in scale_indices) and not (i in exclude_indices and j in exclude_indices):
                scale_entries.add(pair_val - 1)

    except:
        print("Need to supply lennard jones types")
        exit()
# Number of columns in the section
num_columns = 5

# Read all the lines from the input file
with open(input_file) as file:
    lines = file.readlines()

sent = False
ind_counter = 0
Lennard_jones_count = 0

# Define indices that should be scaled

# Read LJ index pairs from the topology file

# Loop through all the lines in the file to find the section
sent = False
ind_counter = 0
for ind, line in enumerate(lines):
    if lookup_str in line:
        first_ind = ind
        sent = True
        continue
    if sent and "FLAG" not in line:
        ind_counter += 1
    elif sent and "FLAG" in line:
        break

# Index at the end of the section
last_ind = first_ind + ind_counter + 1

# Open the output file for writing
with open(output_file, 'w') as file:
    mod_count = 0
    # Loop through all the lines in the file
    for ind, line in enumerate(lines):
        if lookup_str in line:
            file.write(line)
            file.write(lines[ind + 1])
            section = np.array(lines[ind + 2:last_ind])
            new_vals = []

            for sec_line in section:
                for n in sec_line.split():
                    if "\n" in n:
                        n = n[:-1]
                    if lookup_str in ["LENNARD_JONES_ACOEF", "LENNARD_JONES_BCOEF"]:
                        # Apply scaling based on scale_entries, without zero found or float conditions
                        if mod_count in scale_entries:
                            new_vals.append(float(n.lower()) * lamb)
                        else:
                            new_vals.append(float(n.lower()))
                        mod_count += 1
                    elif (lookup_str == "CHARGE" and mod_count < atom_count) or lookup_str == 'DIHEDRAL_FORCE_CONSTANT':                                      
                        new_vals.append(float(n.lower()) * lamb)
                        mod_count += 1
                    else:
                        new_vals.append(float(n.lower()))

            # This is the new line to be saved to the file
            new_line = ''
            for idx, val in enumerate(new_vals):
                # Convert to scientific notation
                if "-" in str(val):
                    sci = ' %.8E' % Decimal(val)
                else:
                    sci = '  %.8E' % Decimal(val)   
                new_line += sci
                if idx % num_columns == 4 or len(new_vals) - 1 == idx:
                    new_line += '\n'
                    file.write(new_line)
                    new_line = ''

            break
        else:
            # Write all the lines on the way to the section
            file.write(line)

    # Write all the lines after breaking
    for l in lines[last_ind:]:
        file.write(l)

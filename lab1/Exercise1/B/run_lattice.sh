#!/bin/bash

# Loop over lattice constants from 3.95 to 4.05 with step size 0.1
for lattice_constant in $(seq 3.95 0.1 4.05); do
    # Generate a temporary LAMMPS input file with the current lattice constant
    cat > lammps_input_temp.in <<EOF
# This section defines the units that you are using. The dimensions of your
# simulation and the periodic boundary conditions.

units metal
dimension 3
boundary p p p

# Define the initial lattice structure and size [angstrom].
atom_style atomic
variable lattice_constant equal ${lattice_constant}
lattice fcc \${lattice_constant}
region box block 0 1 0 1 0 1 # Define box with 1 unit cell per side

# Create the simulation box and populate it with atoms .
create_box 1 box
create_atoms 1 box

# Define atomic mass [g/mol] for the element being simulated .
mass 1 107.8681

# Define the Embedded Atom Potential
pair_style eam
pair_coeff * * Ag_u3.eam

# Define neighbor list settings to improve computational efficiency .
neighbor 0.3 bin
neigh_modify delay 5 every 1

# Allow isotropic box relaxation and atomic position relaxation .
fix 1 all box/relax iso 0.0 vmax 0.001

# Monitor thermodynamic properties with energy output at every step .
thermo 1
thermo_style custom step temp pe lx ly lz press

run 0

# Save the final structure and trajectory .
write_data optimized_structure_${lattice_constant}.data

variable pot_e equal pe
print "\${lattice_constant} \${pot_e}" append energy_volume.dat
EOF

    # Run LAMMPS with the generated input file
    lmp -in lammps_input_temp.in

    # Clean up the temporary input file
    rm lammps_input_temp.in
done

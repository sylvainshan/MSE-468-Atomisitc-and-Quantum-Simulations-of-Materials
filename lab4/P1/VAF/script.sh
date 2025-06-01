#!/bin/bash -f

# Script written to run LAMMPS multiple times on one node.
# The script can loop over 3 different input parameters:
### the time step, specified in list_time_steps
### The temperature, specified in list_temperatures
### The supercell size, specified in list_supercell_size

# There are several ways to set your list in BASH.
### For explicit definition of e.g. your time_steps, you can do
####### list_supercell_size="2 3 4 5 6 7 8 9 10" #######

### You can also use the seq command to create a sequence as in:
####### list_temperatures=`seq 1000 100 3000` #######
### This creates a sequence of values between 1000 and 3000, every 100

list_time_steps="0.001"
list_temperatures="293 800"
list_supercell_size="6"

#################### !IMPORTANT! ####################
#   To run in canonical ensemble for problem 2 comment out line 100 and uncomment line 101
#   Remember to change the 2 to 1 if your virtual machine is running on 1 CPU in line 104
#################### !IMPORTANT! ####################

# This is the executable path for LAMMPS
exec=`which lmp_mpi`

# This is the executable path for the parser
# Change it if you moved the `parser.py` to some other location
parser="python3 /home/max/Desktop/SHARED/LAB4/scripts/parser.py"

# These are constants that the script does not loop over, you need to change them during the exercise
# equilibration and production time in ps, lattice parameter in Angstrom and how often (in ps) everything will be outputted
sim_time="6"
equilibration_time="5"
lattice_parameter="5.37"
stride_time="0.1"

# Start loops
for supercell in $list_supercell_size; do
    for temperature in $list_temperatures; do
        for time_step in $list_time_steps; do
            sim_steps=$(echo "scale=0; $sim_time / $time_step" | bc)
            stride=$(echo "scale=0; $stride_time / $time_step" | bc)
            equilibration_steps=$(echo "scale=0; $equilibration_time / $time_step" | bc)
            echo "Cell size ${supercell}, Temperature ${temperature} K, Timestep $time_step, Equilibration time ${equilibration_steps} steps, Production $sim_steps steps"
            base_name="${supercell}_${temperature}_${time_step}"
            cat > md_${base_name}.in << EOF
# LAMMPS input script for a molecular dynamics simulation of Fe using an EAM potential

# Initialise simulation
clear
# Variables for AgI; after https://aluru.web.engr.illinois.edu/Journals/SSI18.pdf 
# species 1 = (I)
# species 2 = (Ag)
# Conventional cell creation
dimension       3
units           metal
boundary        p p p
atom_style      charge
lattice         custom   1.0     &
        a1      ${lattice_parameter}       0.0     0.0     &
        a2      0.0      ${lattice_parameter}      0.0     &
        a3      0.0      0.0     ${lattice_parameter}      &
        basis   0.00  0.00  0.00              &
        basis   0.50  0.50  0.50              # Lattice for I system
region          box block 0 1 0 1 0 1 units lattice
create_box      2 box # 2 atom-type number (I have two species), box is the name 
create_atoms    1 box # type 1 is I
lattice         custom   1.0     &
        a1      ${lattice_parameter}       0.0     0.0     &
        a2      0.0      ${lattice_parameter}      0.0     &
        a3      0.0      0.0     ${lattice_parameter}      &
        basis   0.25  0.00  0.50              &
        basis   0.75  0.50  0.00              # Lattice for Ag system
create_atoms    2 box # type 2 is Ag
# Masses and charges
mass 1 126.90447      # mass I
mass 2 107.8682       # mass Ag
set type 1 charge -0.3181
set type 2 charge  0.3181
# Ewald sums
kspace_style ewald 1.0e-6  
# Hybrid potential: Morse + Coulomb

pair_style hybrid/overlay coul/long 10.0 morse/smooth/linear 10.0
pair_coeff * * coul/long
#          type type          D      alpha   r0     
pair_coeff   1    2   morse/smooth/linear  0.5500 1.6000 2.6000   
pair_coeff   1    1   morse/smooth/linear  0.1600 0.6840 5.7000

replicate  ${supercell} ${supercell} ${supercell}


# Define simulation parameters
timestep ${time_step}

# Set initial temperature and velocity
velocity all create ${temperature} 87287 loop geom

# Add these lines to check initial velocity distribution
compute vels all property/atom vx vy vz
compute temp_initial all temp
thermo_style custom step temp c_temp_initial
thermo 1
run 0
write_dump all custom initial_velocities_${base_name}.dump id type vx vy vz


# Following would be returned in the output file
thermo ${stride}
thermo_style custom cella cellb cellc
# Run equilibration for the first 5000 timesteps in canonical ensemble
fix 1 all nvt temp ${temperature} ${temperature} $(echo "$time_step * 100.0" | bc)
run ${equilibration_steps}

# Calculate mean square displacement
group I type 1
group Ag type 2
compute mymsdI I msd com yes
compute mymsdAg Ag msd com yes
variable msd_normI equal c_mymsdI[4]
variable msd_normAg equal c_mymsdAg[4]

# Following would be returned in the output file
thermo ${stride}
thermo_style custom time temp pe ke press v_msd_normI v_msd_normAg

# Dumping positions and velocities
dump dp1 all custom ${stride} positions_${base_name}.lammpstrj id type xu yu zu
dump dp2 all custom ${stride} velocities_${base_name}.lammpstrj id type vx vy vz
dump_modify dp2 format line "%d %d %g %g %g" 

# Run production simulation
unfix 1
fix 2 all nve
# fix 2 all nvt temp ${temperature} ${temperature} $(echo "$time_step * 100.0" | bc)
run ${sim_steps}
EOF
            mpirun -np 4 $exec -in md_${base_name}.in > md_${base_name}.out
            $parser md_${base_name}.out positions_${base_name}.lammpstrj velocities_${base_name}.lammpstrj output_${base_name}.json -e E_of_t_${base_name}.dat
        done
    done
done

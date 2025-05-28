#!/usr/bin/bash
# This script calculates the velocity autocorrelation function (VAF) for the specified range of temperatures.
# First run `parser.py` to create the JSON files from LAMMPS output.
# Remeber to change the path to your `vaf.py` script and to modify
# the following variables to your needs.
temps="273 373 473"
supercell="2"
timestep="0.001"
for temp in $temps; do python /your/path/here/vaf.py output_${supercell}_${temp}_${timestep}.json -o vaf_${temp}.dat; done

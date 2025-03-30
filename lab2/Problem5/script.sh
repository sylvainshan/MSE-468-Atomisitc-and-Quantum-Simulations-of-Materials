#!/usr/bin/env bash

# Input data:
LISTA="9.00"      # List of values of lattice parameter to try
LISTECUT="20 30"  # List of plane-wave cutoffs to try
LISTK="2 4"       # List of number of k-points per dimension to try

# Files of interest:
TMP_DIR="./tmp"                 # where temporary data will be stored
PSEUDO_DIR="./pseudopotentials" # where pseudopotentials are stored
OUT_DIR="./simulations"         # where input and output will be
                                # created once the script runs.

PW_LAUNCH='pw.x'         # This is QE executable file.

# Creating directories
if [ ! -d $TMP_DIR ]; then
   mkdir $TMP_DIR
fi

if [ ! -d $OUT_DIR ]; then
   mkdir $OUT_DIR
fi

# Start loops on plane-wave cutoffs, k-point grids, and lattice constants:
for ecut in $LISTECUT; do
   for k in $LISTK; do
      for a in $LISTA; do
      INPUT="$OUT_DIR/CaO.scf.a=$a.ecut=$ecut.k=$k.in"
      OUTPUT="$OUT_DIR/CaO.scf.a=$a.ecut=$ecut.k=$k.out"

      # Create new input file
      cp CaO_primitive.scf.in $INPUT

      sed -i "s/    prefix = .*/    prefix = 'CaO.$a.$ecut.$k'/g" $INPUT
      sed -i "s%    pseudo_dir = .*%    pseudo_dir = '$PSEUDO_DIR'%g" $INPUT
      sed -i "s%    outdir = .*%    outdir = '$TMP_DIR'%g" $INPUT
      sed -i "s/    celldm(1) = .*/    celldm(1) = $a/g" $INPUT
      sed -i "s/    ecutwfc = .*/    ecutwfc = $ecut/g" $INPUT
      sed -i "/K_POINTS/{n;s/.*/    $k $k $k 0 0 0 /}" $INPUT

      # Run PWscf to create new output file:
      echo "Running $PW_LAUNCH < $INPUT > $OUTPUT"
      $PW_LAUNCH < $INPUT > $OUTPUT

      done
   done
done

rm -r $TMP_DIR/*

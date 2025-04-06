#!/usr/bin/env bash

# Input data:
LISTX="0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10"      # Liste des valeurs de déformation (strain x)
LISTECUT="70"                    # Cutoff des ondes planes
LISTK="4"                        # Nombre de k-points par dimension

# Paramètre de réseau à l'équilibre (en Angström)
a0=9.11

# Répertoires
TMP_DIR="./tmp"                 
PSEUDO_DIR="./pseudopotentials" 
OUT_DIR="./tetragonal_strain_c44"   

PW_LAUNCH='mpirun pw.x'  # Exécutable de Quantum ESPRESSO

# Création des répertoires si non existants
mkdir -p $TMP_DIR
mkdir -p $OUT_DIR

# Boucles pour parcourir toutes les combinaisons de paramètres
for ecut in $LISTECUT; do
   for k in $LISTK; do
      for x in $LISTX; do
         INPUT="$OUT_DIR/CaO.scf.x=$x.ecut=$ecut.k=$k.in"
         OUTPUT="$OUT_DIR/CaO.scf.x=$x.ecut=$ecut.k=$k.out"

         # Copier le fichier d'entrée de base
         cp CaO_conventional_c44.scf.in $INPUT

         # Calcul des nouveaux paramètres de réseau
         new_a1=$(echo "$a0 * sqrt(1 + ($x^2)/4)" | bc -l)
         new_a2=$(echo "$a0 * sqrt(1 + ($x^2)/4)" | bc -l)
         new_a3=$(echo "$a0 * (1 + ($x^2)/(4 - ($x^2)))" | bc -l)
         ratio1=$(echo "$new_a2 / $new_a1" | bc -l)
         ratio2=$(echo "$new_a3 / $new_a1" | bc -l)
         angle=$(echo "$x / (1 + ($x^2)/4)" | bc -l)

         # Modifier le fichier d'entrée avec les nouveaux paramètres
         sed -i "s/    prefix = .*/    prefix = 'CaO.$x.$ecut.$k'/g" $INPUT
         sed -i "s%    pseudo_dir = .*%    pseudo_dir = '$PSEUDO_DIR'%g" $INPUT
         sed -i "s%    outdir = .*%    outdir = '$TMP_DIR'%g" $INPUT
         sed -i "s/    celldm(1) = .*/    celldm(1) = $new_a1/g" $INPUT
         sed -i "s/    celldm(2) = .*/    celldm(2) = $ratio1/g" $INPUT
         sed -i "s/    celldm(3) = .*/    celldm(3) = $ratio2/g" $INPUT
         sed -i "s/    celldm(4) = .*/    celldm(4) = $angle/g" $INPUT
         sed -i "s/    ecutwfc = .*/    ecutwfc = $ecut/g" $INPUT
         sed -i "/K_POINTS/{n;s/.*/    $k $k $k 0 0 0 /}" $INPUT

         # Exécuter la simulation
         echo "Running $PW_LAUNCH < $INPUT > $OUTPUT"
         $PW_LAUNCH < $INPUT > $OUTPUT

      done
   done
done

# Nettoyage des fichiers temporaires
rm -r $TMP_DIR/*


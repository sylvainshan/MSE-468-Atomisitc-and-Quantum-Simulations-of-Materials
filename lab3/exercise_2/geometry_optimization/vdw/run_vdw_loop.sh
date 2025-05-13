#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks 18
#SBATCH --cpus-per-task 1
#SBATCH --time=00:10:00
#SBATCH --mem=30000

# Bugfix for oneapi issue
# NOTE: This is only needed on helvetios, uncomment if running there.
# export FI_PROVIDER=verbs

# These are the libraries necessary for our code to run
# NOTE: This is only needed on helvetios, uncomment them if running there.
# module purge
# module load intel/2021.6.0  intel-oneapi-mpi/2021.6.0 quantum-espresso/7.2-mpi

# PW_LAUCHER="srun" # Use this line for helvetios
PW_LAUCHER="mpirun" # Use this for the Virtual Machine

# ---------------
# Initial parameters
# ---------------

PSEUDODIR="../../PP"
OUTODIR='./out'

# ---------------
# Run parameters
# ---------------

THISNAME='MoS2_eos_vdw'

#Note that in this case you can optimize A and C separatly 
#at first put only the experimental value of C n CLIST and optimize A
#then put the optimized value of A in ALIST and optimize C
#you can optimize them together if you want but be careful with the number of caluclations
#if you put 10 values both in ALIST and 10 in CLIST you will end up with 100 calcs.
#The values in the lists below will be used for the parameter 'celldm' in the QuantumESPRESSO input file. Check on the QE documentation what is the correct unit for this parameter.

AEXP="5.9655"
CEXP="23.243631348665723"

startA=$(echo "$AEXP-0.1" | bc -l)
endA=$(echo "$AEXP+0.1" | bc -l)

startC=$(echo "$CEXP-0.2" | bc -l)
endC=$(echo "$CEXP+0.2" | bc -l)

#Note, please uncomment one of the two lines below to reflex what you are doing

# OPTIMIZE="optA"
OPTIMIZE="optC"

if [ "$OPTIMIZE" = "optA" ]; then
  ALIST=$(seq $startA 0.01 $endA)
  CLIST=$CEXP
elif [ "$OPTIMIZE" = "optC" ]; then
  ALIST=$AEXP
  CLIST=$(seq $startC 0.01 $endC)
else
  echo "Erreur: OPTIMIZE doit être 'optA' ou 'optC'"
  exit 1
fi

# ---------------
# Running the calculations
# ---------------

OUT_SUMMARY_NAMEac=$THISNAME"_ac_"$OPTIMIZE".dat"
OUT_SUMMARY_NAMEacoa=$THISNAME"_acoa_"$OPTIMIZE".dat"

echo "# a    c    total energy"   >> $OUT_SUMMARY_NAMEac
echo "# a   c/a   total energy"   >> $OUT_SUMMARY_NAMEacoa

for A in $ALIST
do
for C in $CLIST
do

CoA=`echo "scale=8;$C/$A" | bc`

inputname=$THISNAME\_a=$A\_c=$C\_scf.in
outputname=$THISNAME\_a=$A\_c=$C\_scf.out

echo "Running $inputname"

cat <<EOF > $inputname 
&control
  calculation  = 'relax'
  restart_mode = 'from_scratch'
  prefix       = '$THISNAME'
  pseudo_dir   = '$PSEUDODIR'
  outdir       = '$OUTODIR'
 /
&system
  ibrav=4, 
  celldm(1)=$A
  celldm(3)=$CoA 
  nat=  6, 
  ntyp= 2,
  ecutwfc = 40
  ecutrho = 600
  input_dft = ... complete this line with instructions in the handouts.
/
&ELECTRONS
  conv_thr =   1.0000000000d-8,
  electron_maxstep = 200,
  mixing_beta =   3.0000000000d-01,
  mixing_mode = 'plain',
  startingwfc = 'atomic+random',
/
&IONS
/
K_POINTS automatic
4 4 2 0 0 0

ATOMIC_SPECIES
Mo  95.940000 Mo.pbe-spn-rrkjus-tested-pslib025.UPF
S   32.066000 S.pbe-n-rrkjus-tested-pslib025.UPF

ATOMIC_POSITIONS crystal
Mo            0.3333333333        0.6666666667        0.0000000000
S             0.6666666667        0.3333333333       -0.1283498612
S             0.6666666667        0.3333333333        0.1283498612
Mo            0.6666666667        0.3333333333        0.5000000000
S             0.3333333333        0.6666666667        0.3716501388
S             0.3333333333        0.6666666667        0.6283498612
EOF

$PW_LAUCHER pw.x < $inputname > $outputname 2> /dev/null

EN=`cat $outputname | grep -e ! | egrep -o "([+-])?[0-9]+(\.[0-9]+)?" | tail -1`

echo "$A $C $EN"   >> $OUT_SUMMARY_NAMEac
echo "$A $CoA $EN" >> $OUT_SUMMARY_NAMEacoa

echo "Ok - Energy $EN"
done
done


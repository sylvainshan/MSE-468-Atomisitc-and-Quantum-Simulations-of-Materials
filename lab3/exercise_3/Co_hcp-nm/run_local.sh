#!/usr/bin/env bash

if [[ -z $1 ]]; then

  echo '❗️Please specify an input file!'

elif [[ ! -e $1 ]]; then

  echo '❗️Specified file does not exist'

else

  if [[ -n $(grep '< TO COMPLETE >' $1) ]]; then
    echo '🤨 Did you fill in all the `< TO COMPLETE >` inputs?'
  else
    # Figure out the executable to run
    if [[ -n $(grep -i '&control' $1) ]]; then executable=pw.x
    elif [[ -n $(grep -i '&dos' $1) ]]; then executable=dos.x
    elif [[ -n $(grep -i '&projwfc' $1) ]]; then executable=projwfc.x
    elif [[ -n $(grep -i '&bands' $1) ]]; then executable=bands.x
    fi

    # Bombs away!
    mpirun $executable -in $1 > ${1%.in}.out 2> /dev/null

    EN=` grep ! scf.out | egrep -o "([+-])?[0-9]+(\.[0-9]+)?" `
    if [[ -n $EN ]]; then echo "Energy = ${EN} Ry"; fi
  fi
fi

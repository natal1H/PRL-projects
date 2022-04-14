#!/bin/bash

if [ $# -eq 1 ]; then
    proc=2

    #preklad cpp zdrojaku
    #mpic++ --prefix /usr/local/share/OpenMPI -o pro pro.cpp
    mpic++ --prefix /usr/share/OpenMPI -o pro pro.cpp

    #spusteni
    #mpirun -oversubscribe --prefix /usr/local/share/OpenMPI -np $proc pro
    mpirun --prefix /usr/share/OpenMPI -np $proc pro $1

    #uklid
    rm -f pro
else
    echo "Nespravne parametre. Spustenie: ./test.sh RETAZEC"
fi;


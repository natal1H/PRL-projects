#!/bin/bash

numbers=8
proc=19

#preklad cpp zdrojaku
mpic++ --prefix /usr/local/share/OpenMPI -o oems oems.cpp

#vyrobeni souboru s random cisly
dd if=/dev/random bs=1 count=$numbers of=numbers

#spusteni
mpirun -oversubscribe --prefix /usr/local/share/OpenMPI -np $proc oems

#uklid
rm -f oems numbers
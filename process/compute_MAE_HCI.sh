#!/bin/bash

nstate=70
nmax=2
ntot=2
eps2=0
perturb=0
molecule="acetonitrile"
ref_file="ref_${molecule}_AVCI"
rm varstates.dat -f
rm ptstates.dat -f
rm maxerrors.dat -f
rm rmserrors.dat -f

rm test_energies.dat -f
rm ref_energies.dat -f
rm temp.dat -f
rm errors.dat -f
rm errors_sq.dat -f
for state in `seq 1 ${nstate}`; do
    awk -v VAR1=$state '(FNR==VAR1){printf "%0.3f ", $1;printf "\n";}' cipsi_0.00022.dat >> test_energies.dat
    awk -v VAR1=$state '(FNR==(VAR1)){printf "%0.3f ", $1;printf "\n";}' ${ref_file}.dat >> ref_energies.dat
done
paste test_energies.dat ref_energies.dat > temp.dat
awk '($1>=$2){print $1-$2}($1<$2){print $2-$1}' temp.dat > errors.dat
awk '{print ($1-$2)^2}' temp.dat > errors_sq.dat
awk 'BEGIN{a=0}{if ($1>0+a) a=$1} END{print a}' errors.dat >> maxerrors.dat
awk '{ sum += $1 } END { if (NR > 0) print sqrt(sum / NR)}' errors_sq.dat >> rmserrors.dat

paste maxerrors.dat rmserrors.dat >> MAE_HCI_CIPSI_0.00022.dat

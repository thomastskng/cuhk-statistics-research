#!/bin/sh
#trial.sh
#Torque script to run Matlab program

#Torque directives
#PBS -N QRGB-LB-n200m250t50c20NormalCase-Naive
#PBS -l nodes=1,walltime=12:00:00,mem=250mb
#PBS -m abe
#PBS -V

cd $PBS_O_WORKDIR/sim$PBS_ARRAYID

#for i in ${PBS_ARRAYID[@]}; do
#	echo "${i}"
#done
unique_moneyness=(0.80 0.85 0.90 0.95 0.96 0.97 0.98 0.99 1.00)
moneyness=($(for i in ${unique_moneyness[@]}; do for ((n=1; n<=10; n++)) do echo -n "$i ";done;done;))

#set output and error directories
#PBS -o $PBS_O_WORKDIR/sim$PBS_ARRAYID
#PBS -e $PBS_O_WORKDIR/sim$PBS_ARRAYID

#define parameter lambda
export LAMBDA=10

#Command to execute Matlab code
# matlab -nosplash -nodisplay -nodesktop -r emleoption_N1000K4I500PBSMJ $PBS_ARRAYID  > matoutfile
source /opt/share/vni/imsl/fnl600/lnxin100e64/bin/fnlsetup.sh
#Command below is to execute Matlab code for Job Array (Example 4) so that each part writes own output
matlab -nosplash -nodisplay -nodesktop -r "stat(${moneyness[$PBS_ARRAYID]},$PBS_ARRAYID)" $PBS_ARRAYID  > matoutfile.$PBS_ARRAYID

#!/bin/bash -login
#PBS -l walltime=35:00:00
#PBS -l nodes=1:ppn=20
#PBS -l mem=100gb
#PBS -A hirn
#PBS -j oe
#PBS -N SuccessGrid_k5_n20_rhoOne

cd ${PBS_O_WORKDIR}

### Load MATLAB
#PBS -W x=gres:MATLAB

### Change MATLAB version
module swap MATLAB/R2014a MATLAB/R2017b

# call your executable
matlab-mt -nodisplay -r "SuccessGrid_k5_n20_rhoOne"
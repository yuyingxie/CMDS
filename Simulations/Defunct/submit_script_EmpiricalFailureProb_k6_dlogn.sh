#!/bin/bash -login
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=10
#PBS -l mem=64gb
#PBS -A hirn
#PBS -j oe
#PBS -N EmpiricalFailureProb_k6_dlogn

cd ${PBS_O_WORKDIR}

### Load MATLAB
#PBS -W x=gres:MATLAB

### Change MATLAB version
module swap MATLAB/R2014a MATLAB/R2017b

# call your executable
matlab-mt -nodisplay -r "EmpiricalFailureProb_k6_dlogn"
#!/bin/bash -login
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --job-name k3_n20_rhoLarge_DenoiseEigs

cd ~/MDScode

### Change MATLAB version
module load MATLAB/2017b

# call your executable
srun matlab-mt -nodisplay -r "SuccessGrid_k3_n20_rhoLarge_DenoiseEigs"
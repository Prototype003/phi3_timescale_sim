#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Kai Xi,  Oct 2014
#              help@massive.org.au
# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'
# $1: line counter
# Need to use variables OUTSIDE of this script, #SBATCH doesn't support variables: https://help.rc.ufl.edu/doc/Using_Variables_in_SLURM_Jobs

#SBATCH --job-name=compute

# To set a project account for credit charging, 
#SBATCH --account=ot95

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
# SBATCH --exclusive
# SBATCH --cpus-per-task=16

# Memory usage (MB)
#SBATCH --mem-per-cpu=500

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-03:00:00

# SBATCH --qos=shortq
# SBATCH --partition=short,comp

# To receive an email when job completes or fails
# SBATCH --mail-user=aleu6@student.monash.edu
# SBATCH --mail-type=END
# SBATCH --mail-type=FAIL

# Set the file for output (stdout)
#SBATCH --output=compute.out

# Set the file for error log (stderr)
#SBATCH --error=compute.err

# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name

# Job script

# Compute TPMs
# module load matlab/r2019a
# pushd ../
# matlab -nodisplay -nodesktop -r "main_tpms_tau_multipleSystems; exit"
# echo "tpms built, now tarring"
# ./tpm_tar.bash 3chMotifsNLThresh0_9Lag9-11_nSamples200000_nRuns10_medianSplit_tauSearch_nCh2
# echo "tpms tarred, now renaming tar"
# cd tpms/3chMotifsNLThresh0_9Lag9-11_nSamples200000_nRuns10_medianSplit_tauSearch_nCh2/
# mv 3chMotifsNLThresh0_9Lag9-11_nSamples200000_nRuns10_medianSplit_tauSearch_nCh2.tar tpms.tar
# echo "tar renamed"
# cd ../
# popd

# Compute phis
module load python/3.6.2
source ../pyphi_environment/bin/activate
time python phi_compute.py split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_medianSplit_tauSearch_binAverage
deactivate

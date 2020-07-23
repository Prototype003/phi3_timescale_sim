#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Kai Xi,  Oct 2014
#              help@massive.org.au
# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'
# $1: line counter
# Need to use variables OUTSIDE of this script, #SBATCH doesn't support variables: https://help.rc.ufl.edu/doc/Using_Variables_in_SLURM_Jobs

#SBATCH --job-name=loopAll

# To set a project account for credit charging, 
#SBATCH --account=ot95

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
# SBATCH --exclusive
#SBATCH --cpus-per-task=14

# Memory usage (MB)
#SBATCH --mem-per-cpu=1000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-10:00:00
# SBATCH --qos=shortq
# SBATCH --partition=short,comp

# To receive an email when job completes or fails
# SBATCH --mail-user=aleu6@student.monash.edu
# SBATCH --mail-type=END
# SBATCH --mail-type=FAIL

# Set the file for output (stdout)
#SBATCH --output=loopAll.out

# Set the file for error log (stderr)
#SBATCH --error=loopAll.err

# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name

# Job script
module load matlab/r2019b
time matlab -nodisplay -nodesktop -r "main_tpms_thresh_data_stateCountMatch; exit"

module load python/3.6.2
source pyphi_environment/bin/activate
cd phi_3/
time python phi_compute.py split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_threshSplit_binAverage_100perState

deactivate

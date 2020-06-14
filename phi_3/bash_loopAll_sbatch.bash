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
#SBATCH --cpus-per-task=1
# Memory usage (MB)
#SBATCH --mem-per-cpu=2000
# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-06:00:00
# SBATCH --qos=shortq
# SBATCH --partition=short,comp
# To receive an email when job completes or fails
# SBATCH --mail-user=aleu6@student.monash.edu
# SBATCH --mail-type=END
# SBATCH --mail-type=FAIL
# Set the file for output (stdout)
#SBATCH --output=logs/loopAll.out
# Set the file for error log (stderr)
#SBATCH --error=logs/loopAll.err
# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name
# Job script
module load python/3.6.2
source ../pyphi_environment/bin/activate
lines=$(wc -l < array_commands) # Total number of jobs which need to be computed
line_increment=1
# Loop through parameter lines
for (( line=1; line<=$lines; line=$line+$line_increment )); do
	echo "line $line"
	command=$(sed -n "$((${line}))p" "array_commands")
	if [[ -n "$command" ]]; then
		time eval $command
	fi
	echo "done"
done
deactivate

#!/bin/bash
# Loop through everything on own machine, without MASSIVE/Slurm

# Activate virtual environment
source ../pyphi_environment/bin/activate

lines=$(wc -l < array_commands) # Total number of jobs which need to be computed
line_increment=1
# Loop through parameter lines
for (( line=1; line<=$lines; line=$line+$line_increment )); do
	echo "line $line"
	command=$(sed -n "$((${line}))p" "array_commands")
	if [[ -n "$command" ]]; then
		time eval $command > log_local
	fi
	echo "done"
done
deactivate

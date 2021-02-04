#!/bin/bash
nChannels=2
flies=1 # 13
conditions=(1)
#set_ids=($(seq 1 1 1365)) # ($(seq first step last)) # (1036)
set_ids=($(seq 1 1 105))
#set_ids=(1)
#taus=(1 2 3 4 5 10 20 30 40 50 75 100 125 150 175 200 225 250) #(4)
taus=($(seq 1 1 11))
trials=1
global_tpm=1
tau_bin=1
sample_offset=1
bin_method="diff"
> array_commands
for (( fly=1; fly<=$flies; fly++ )); do
	for condition in "${conditions[@]}"; do
		for set in "${set_ids[@]}"; do
			for tau in "${taus[@]}"; do
				for (( trial=1; trial<=$trials; trial++ )); do
					
					echo "python3 phi_compute.py $nChannels $fly $condition $set $tau $trial $global_tpm $tau_bin $sample_offset $bin_method" >> array_commands
					
				done
			done
		done
	done
done

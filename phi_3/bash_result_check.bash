#!/bin/bash
nChannels=4
flies=13 # 13
conditions=(1 2)
set_ids=($(seq 1 1 1365)) #(1036)
taus=(4) #(1 2 3 4 8 12 16 24 32 48 64 128 256 512 4500 9000) #(4 8 16)
trials=8
global_tpm=0
tau_bin=0
sample_offset=0
prefix="split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials_nComponents4_ic2channels_nChannels${nChannels}_globalTPM${global_tpm}_"
suffix=".mat"
if [ $global_tpm -eq 0 ]; then
	tau_string="tau"
	offset_string=""
else
	tau_string="tauBin"
	offset_string="binOffset${sample_offset}"
fi
> array_commands
for (( fly=1; fly<=$flies; fly++ )); do
	printf -v fly_padded "%02d" $fly
	for condition in "${conditions[@]}"; do
		for set in "${set_ids[@]}"; do
			printf -v set_padded "%04d" $set
			for tau in "${taus[@]}"; do
				for (( trial=1; trial<=$trials; trial++ )); do
					
					id="f${fly_padded}c${condition}${tau_string}${tau}${offset_string}s${set_padded}t${trial}"
					
					results_file="results_split_ic2channels/${prefix}${id}${suffix}"
					
					# Check if results file exists, if not, output command to recompute
					if [ ! -e $results_file ]; then
						echo "${results_file}"
						echo "python3 phi_compute.py $nChannels $fly $condition $set $tau $trial $global_tpm $tau_bin $sample_offset" >> array_commands
					fi
					
				done
			done
		done
	done
done

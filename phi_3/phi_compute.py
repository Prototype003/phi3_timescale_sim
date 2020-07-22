
import os
import sys
import tarfile
import multiprocessing as mp
sys.path.append('../')
import numpy as np
import scipy.signal as sp_signal
import scipy.stats as sp_stats
import pyphi
import itertools
from fly_phi import *

def tar_compute(tpm_file, tpm_dir, tpm_type, results_directory, queue):
	tpm_archive = tarfile.open(tpm_dir+tpm_type+'.tar', mode='r')
	
	# Extract files
	tpm_archive.extract(tpm_file, path=tpm_dir)
	
	# Load TPM
	loaded = load_mat(tpm_dir + tpm_file.name)
	tpm = loaded['tpm']
	state_counters = loaded['state_counters']
	nValues = loaded['nValues'][0][0]
	tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm) # assumes loaded TPMs are state-by-state
	
	# Compute phi
	out_file = phi_compute(tpm_formatted, state_counters, nValues, results_directory, tpm_file.name, tpm_type)
	
	# Add output to queue
	queue.put(out_file)
	print(out_file)
	
	# Delete extracted file
	os.remove(tpm_dir+tpm_file.name)


def phi_compute(tpm, state_counters, nValues, out_dir, out_file, tpm_type):
	# Assumes tpm is state-by-node
	# Intputs:
	#	tpm = state-by-node TPM
	#	state_counters = vector of state counts
	#	nValues = number of states each element can take
	#	out_dir = string; directory to output results file
	#	out_name = string; name of output files
	
	# Build pyphi network
	network = pyphi.Network(tpm)

	# Determine number of system elements
	nChannels = np.shape(tpm)[1]
	
	# Determine number of system states
	n_states = nValues ** nChannels

	# Results structure
	phi = dict()

	# Initialise results storage structures
	phi_value = 0; # Doesn't need initialisation
	mips = np.empty((n_states), dtype=tuple)
	big_mips = np.empty((n_states), dtype=object)
	state_phis = np.zeros((n_states))

	# sys.exit()

	# Calculate all possible phi values (number of phi values is limited by the number of possible states)
	for state_index in range(0, n_states):
		#print('State ' + str(state_index))
		# Figure out the state
		state = pyphi.convert.loli_index2state(state_index, nChannels)
		
		# As the network is already limited to the channel set, the subsystem would have the same nodes as the full network
		subsystem = pyphi.Subsystem(network, state, network.node_indices)
		
		#sys.exit()
		
		# Compute phi values for all partitions
		big_mip = pyphi.compute.big_mip(subsystem)
		
		# Store phi and associated MIP
		state_phis[state_index] = big_mip.phi
		mips[state_index] = big_mip.cut
		
		# MATLAB friendly storage format (python saves json as nested dict)
		big_mip_json = big_mip.to_json()
		# Sort each constellation by their mechanisms
		big_mip_json['partitioned_constellation'] = sort_constellation_by_mechanism(big_mip_json['partitioned_constellation'])
		big_mip_json['unpartitioned_constellation'] = sort_constellation_by_mechanism(big_mip_json['unpartitioned_constellation'])
		# Store big_mip
		big_mips[state_index] = big_mip_json
		
		print('State ' + str(state_index) + ' Phi=' + str(big_mip.phi))

	phi_total = 0
	for state_index in range(0, n_states):
		if state_counters[state_index] > 0:
			# Add phi to total, weighted by the number of times the state occurred
			phi_total += state_phis[state_index] * state_counters[state_index]
	phi_value = phi_total / np.sum(state_counters)

	phi['phi'] = phi_value
	phi['state_counters'] = state_counters
	phi['big_mips'] = big_mips
	phi['state_phis'] = state_phis
	phi['tpm'] = tpm
	phi['mips'] = mips

	# Save
	save_mat(out_dir+out_file, {'phi': phi})
	print('saved ' + out_dir + out_file)
	
	return out_file
	

def tar_append(out_dir, queue):
	# Open tar for appending
	archive = tarfile.open(out_dir+'phis.tar', mode='a')
	print('tar open')
	
	# Wait for other processes to finish, and append to tar
	while 1:
		out_file = queue.get()
		print(out_file)
		if out_file == 'finish': # Check if finished - if so, then exit
			archive.close()
			break
		archive.add(out_dir+out_file)
		print('added to tar:'+out_file)
		
		# Remove associated mat file
		os.remove(out_dir+out_file)
		print(out_dir+out_file)
	

# Setup ############################################################################

# Remove PyPhi parallelisation (as we have our own parallelisation)
pyphi.config.PARALLEL_CONCEPT_EVALUATION = False
pyphi.config.PARALLEL_CUT_EVALUATION = False

# Parameters for loading data and calculating
tpm_type = sys.argv[1]

# TPM locations
tpm_dir = "../tpms/" + tpm_type + "/"

# Output location
results_directory = "results/split/" + tpm_type + "/"
if not os.path.exists(results_directory):
	os.makedirs(results_directory)

# Loop through TPMs ############################################################################

# Check for tar archive 
if tarfile.is_tarfile(tpm_dir+tpm_type+'.tar'):
	tpm_archive = tarfile.open(tpm_dir+tpm_type+'.tar', mode='r')
	tpms = tpm_archive.getmembers()
	
	# Create pool
	manager = mp.Manager()
	queue = manager.Queue()
	pool = mp.Pool(mp.cpu_count())
	
	# Create process for appending to tar
	writer = pool.apply_async(tar_append, (results_directory, queue))
	
	# Iterate through TPM files in the tar archive
	jobs = []
	for tpm_file in tpms:
		job = pool.apply_async(tar_compute, (tpm_file, tpm_dir, tpm_type, results_directory, queue))
		jobs.append(job)
	
	# Force wait until all jobs finish
	for job in jobs:
		job.get()
	
	queue.put('finish') # Signal to writer to finish up
	
	# Close pool
	pool.close()
	pool.join()
	
else:
	# Iterate through separate individual TPM files
	for tpm_file in os.listdir(tpm_dir):
		print(tpm_file)
		
		if tpm_file != "params.mat":
			# Check if phi result already exists
			if not os.path.isfile(results_directory+tpm_file):
				# Load TPM
				loaded = load_mat(tpm_dir + tpm_file)
				tpm = loaded['tpm']
				state_counters = loaded['state_counters']
				nValues = loaded['nValues'][0][0]
				tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm) # assumes loaded TPMs are state-by-state
				
				phi_compute(tpm_formatted, state_counters, nValues, results_directory, tpm_file, tpm_type)

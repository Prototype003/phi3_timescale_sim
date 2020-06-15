
import os
import sys
sys.path.append('../')
import numpy as np
import scipy.signal as sp_signal
import scipy.stats as sp_stats
import pyphi
import itertools
from fly_phi import *

# Calculate for only 1 fly, 1 set size, etc (see input parameters below) - loop using bash script

# Setup ############################################################################

# Parameters for loading data and calculating
nChannels = int(sys.argv[1]) # Number of channels across which to compute phi
fly = int(sys.argv[2]) # 1-indexed
condition = int(sys.argv[3]) # 1-indexed
set = int(sys.argv[4]) # 1-indexed
tau = int(sys.argv[5]) # Actual tau in terms of how many samples to lag across (in ms, for 1000Hz sampling rate)
trial = int(sys.argv[6]) # 1-indexed
global_tpm = int(sys.argv[7]) # 0=8 trials (2250 samples); 1=1 trial (18000 samples)
tau_bin = int(sys.argv[8]) # 0=don't average across tau samples, use stepsize tau; 1=average across tau samples, afterwards use stepsize 1
sample_offsets = int(sys.argv[9]) # Compute TPM across sample offsets before binning (for when tau_bin == 1)
#binarise_method = sys.argv[10] # Method for binarising - 'median' or 'diff'
binarise_method = 'median'

# Fly data location
data_directory = "../sim_data/"
data_file_prefix = "unidir_no_inst_nSamples18000_nRuns10"
data_file = data_file_prefix + ".mat"

# Output location
results_directory = "results_split_nRuns10/"
if not os.path.exists(results_directory):
	os.makedirs(results_directory)

# tau string for results file
if tau_bin == 1:
	tau_type = "tauBin"
	tau_string = tau_type + str(tau) + "binOffset" + str(sample_offsets)
else:
	tau_type = "tau"
	tau_string = tau_type + str(tau)

# Results file
if binarise_method == 'diff':
	results_file_suffix = "_nChannels" + str(nChannels) + "_" + binarise_method + "_globalTPM" + str(global_tpm) + "_f" + "{0:0>2}".format(fly) + "c" + str(condition) + tau_string + "s" + "{0:0>4}".format(set) + "t" + str(trial)
else:
	results_file_suffix = "_nChannels" + str(nChannels) + "_globalTPM" + str(global_tpm) + "_f" + "{0:0>2}".format(fly) + "c" + str(condition) + tau_string + "s" + "{0:0>4}".format(set) + "t" + str(trial)
results_file = data_file_prefix + results_file_suffix + ".mat"

# Load data ############################################################################

loaded_data = load_mat(data_directory + data_file)
fly_data = loaded_data['data']
print("Fly data loaded")

# Filter for data from channel set of specific fly
channel_combinations = list(itertools.combinations(np.arange(fly_data.shape[1]), nChannels)) # All channel combinations
channel_set = np.asarray(channel_combinations[set-1])

fly_data = fly_data[:, :, :, :, None] # Format data to match fly data (5D - samples x channels x trials x flies x conditions)
fly_data = fly_data[:, :, :, fly-1, condition-1] # From here, only holds data for specified fly, condition
fly_data = fly_data[:, channel_set, :] # Data for channel set (can't do with fly,condition because dimension order changes for some reason)

# Get trial
if global_tpm == 1:
	# Flatten trials if using global TPM
	flattened_data = fly_data[:, :, 0]
	for trial in range(1, fly_data.shape[2]):
		flattened_data = np.concatenate((flattened_data, fly_data[:, :, trial]), axis=0)
	fly_data = flattened_data
else:
	# Get specific trial
	fly_data = fly_data[:, :, trial-1] # From here, only holds data for specified trial

print("Specific data obtained")

# Preprocess ############################################################################

tpm_built = 0
if tau_bin == 1:
	if sample_offsets == 1:
		n_values = 2
		tpm, transition_counters = build_tpm_bin_offsets(fly_data, n_values, tau, binarise_method)
		tpm_formatted = tpm
		tpm_built = 1
	else:
		# Downsample by averaging in bins of length tau
		fly_data = tau_resample(fly_data[:, :, None, None, None], tau)
		fly_data = fly_data[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions
		tau_step = 1
else:
	tau_step = tau

if tpm_built == 0:
	# Binarise by median split
	fly_data, n_values, medians = binarise_trial_median(fly_data[:, :, None, None, None])
	fly_data = fly_data[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions

	# Build TPM
	tpm = build_tpm(fly_data[:, :, None], tau_step, n_values)
	tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm)

print("TPM built")

# if sample_offsets == 1:
	# n_values = 2
	# tpm, transition_counters = build_tpm_bin_offsets(fly_data, n_values, tau)
	# tpm_formatted = tpm
# else:
	# if tau_bin == 1:
		# # Downsample by averaging in bins of length tau
		# fly_data = tau_resample(fly_data[:, :, None, None, None], tau)
		# fly_data = fly_data[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions
		# tau_step = 1
	# else:
		# tau_step = tau
	
	# # Binarise by median split
	# fly_data, n_values, medians = binarise_trial_median(fly_data[:, :, None, None, None])
	# fly_data = fly_data[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions

	# # Build TPM
	# tpm = build_tpm(fly_data[:, :, None], tau_step, n_values)
	# tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm)
	# #tpm = build_tpm(np.flip(fly_data[:, :, None], 0), tau_step, n_values)
	# #tpm, tmp = build_tpm_sbn_normalise(fly_data[:, :, None], tau_step, n_values, 9000)
# print("TPM built")

# Build the network and subsystem
# We are assuming full connection
network = pyphi.Network(tpm_formatted)
print("Network built")

#########################################################################################
# Remember that the data is in the form a matrix
# Matrix dimensions: sample(2250) x channel(15)

# Determine number of system states
n_states = n_values ** len(channel_set)

# Results structure
phi = dict()
phi['nChannels'] = nChannels
phi['channel_set'] = channel_set + 1 # Save 1-indexed values
phi['tau'] = tau

# Initialise results storage structures
phi_value = 0; # Doesn't need initialisation
mips = np.empty((n_states), dtype=tuple)
big_mips = np.empty((n_states), dtype=object)
state_counters = np.zeros((n_states))
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

# We average across the phi values previously computed for the channel set
# When averaging across samples, we weight by the number of times each state occurred
if sample_offsets == 0:
	for sample_counter in range(0, fly_data.shape[0]):
		sample = fly_data[sample_counter, :]
		
		# Determine the state
		state = pyphi.convert.state2loli_index(tuple(sample))
		
		# Add to state_counter
		state_counters[state] += 1
else:
	state_counters = transition_counters

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

# Save ###########################################################################

save_mat(results_directory+results_file, {'phi': phi})
print('saved ' + results_directory + results_file)

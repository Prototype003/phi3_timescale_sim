
import pyphi as pyphi
import numpy as np
import scipy.io as sio
#import matplotlib.pyplot as plt
#import sklearn.preprocessing as skp
#import copy
#import itertools

print("inside")

# Load functions
def load_mat(data_location):
	# data_location should give the directory location and file name (i.e. the path to the file)
	import scipy.io as sio
	return sio.loadmat(data_location)

def save_mat(data_directory, dictionary):
	import scipy.io as sio
	sio.savemat(data_directory, dictionary, do_compression=True)

def tau_resample(fly_data, tau):
	"""
	Averages across every tau samples
	Any remaining samples are dropped
	Source - https://stackoverflow.com/questions/30379311/fast-way-to-take-average-of-every-n-rows-in-a-npy-array
	Inputs:
		fly_data = data matrix with dimensions (samples x channels x trials x flies x conditions)
		tau = number of samples to average across
	Outputs:
		fly_data_resampled = matrix with same dimensions as fly_data, except for samples, which has size (samples / tau)
	"""
	
	cum = np.cumsum(fly_data, 0) # Sum cumulatively across samples
	fly_data_resampled = cum[tau-1::tau, :, :, :, :] / float(tau) # Take every tau'th cumulative sum, and divide by tau
	fly_data_resampled = fly_data_resampled[1:, :, :, :, :] - fly_data_resampled[:-1, :, :, :, :] # Subtract previous 'cumulative average'
	
	return fly_data_resampled

def binarise(value, threshold):
	"""
	Returns 1 if value > threshold
	Returns 0 if value <= threshold
	Inputs:
		value = some number
		threshold = some number
	Outputs:
		binarised = 1 or 0
	"""
	
	if value > threshold:
		binarised = 1
	else:
		binarised = 0
	
	return binarised

def binarise_global_median(fly_data):
	"""
	Finds the median for each channel, across all epoch-trials
	Binarises samples based on this median: 1 if greater than median, 0 otherwise
	Inputs:
		fly_data = data matrix with dimensions (samples x channels x epoch-trials x flies x conditions)
	Outputs:
		fly_data_binarised = matrix with same dimensions as fly_data, with binarised values
		channel_medians = matrix with medians used for binarisation
	"""
	
	lengths = np.shape(fly_data)
	
	# Reorder dimensions to prepare for reshape
	fly_data_reordered = np.swapaxes(fly_data, 1, 2)
	# Reformat - collapse on epoch-trials so that first dimension is samples x epoch-trials
	# Option 'F' means the first index changes fastest, second index changes sexond fastest, etc.
	fly_data_globalised = np.reshape(fly_data_reordered, (lengths[0]*lengths[2], lengths[1], lengths[3], lengths[4]), 'F')
	
	# Get median per channel, fly, condition
	channel_medians = np.median(fly_data_globalised, axis=0)
	
	# Binarise based on median:
	# Greater than median = 1
	# Less than or equal to median = 0
	fly_data_binarised = np.ndarray(np.shape(fly_data_globalised))
	binarise_vector = np.vectorize(binarise) # We will binarise channel (vector of samples) by channel, but this can probably be done cleaner through Generalized Universal Functions
	
	for condition in range(0, np.shape(fly_data_globalised)[3]):
		for fly in range(0, np.shape(fly_data_globalised)[2]):
			for channel in range(0, np.shape(fly_data_globalised)[1]):
				fly_data_binarised[:, channel, fly, condition] = binarise_vector(fly_data_globalised[:, channel, fly, condition], channel_medians[channel, fly, condition])
	
	# Return to original shape
	fly_data_binarised = np.split(fly_data_binarised, lengths[2], 0) # This splits the globalised samples into trials (re-adding the trials dimension as the first dimension)
	fly_data_binarised = np.transpose(fly_data_binarised, (1, 2, 0, 3, 4))
	
	return fly_data_binarised, 2, channel_medians

def binarise_trial_median(fly_data):
	"""
	Finds the median for each channel, at each epoch-trial
	Binarises samples based on this median: 1 if greater than median, 0 otherwise
	Inputs:
		fly_data = data matrix with dimensions (samples x channels x epoch-trials x flies x conditions)
	Outputs:
		fly_data_binarised = matrix with same dimensions as fly_data, with binarised values
		channel_medians = matrix holding medians used for binarisation
	"""
	
	# Get channel medians across samples, per channel, trial, fly, and condition
	channel_medians = np.median(fly_data, axis=0)
	
	# Binarise based on median:
	# Greater than median = 1
	# Less than or equal to median = 0
	fly_data_binarised = np.ndarray(np.shape(fly_data))
	binarise_vector = np.vectorize(binarise) # We will binarise channel (vector of samples) by channel, but this can probably be done cleaner through Generalized Universal Functions
	
	for condition in range(0, np.shape(fly_data)[4]):
		for fly in range(0, np.shape(fly_data)[3]):
			for trial in range(0, np.shape(fly_data)[2]):
				for channel in range(0, np.shape(fly_data)[1]):
					fly_data_binarised[:, channel, trial, fly, condition] = binarise_vector(fly_data[:, channel, trial, fly, condition], channel_medians[channel, trial, fly, condition])
	
	return fly_data_binarised, 2, channel_medians

def binarise_trial_diff(fly_data):
	"""
	Binarises time series baesd on gradient between adjacent samples: 1 if increasing, 0 otherwise
	Inputs:
		fly_data = data matrix with dimensions (samples x channels x epoch-trials x flies x conditions)
	Outputs:
		fly_data_binarised = matrix with same dimensions as fly_data, with binarised values
			Samples dimension will be 1 element smaller
		channel_gradients = matrix holding gradients used for binarisation
	"""
	
	# Get gradients across samples, per channel, trial, fly, and condition
	channel_gradients = np.diff(fly_data, n=1, axis=0)
	
	# Binarise baesd on gradients:
	# Increasing = 1
	# Decreasing or unchanging = 0
	fly_data_binarised = (channel_gradients > 0).astype(int)
	
	return fly_data_binarised, 2, channel_gradients
	

def build_tpm(fly_data, tau, n_values):
	"""
	Builds a TPM for one fly and one condition, holding the probabilities that each combination of
	channel-states will (i.e. each system state) will transition into another set of channel-states
	
	Inputs:
		fly_data = matrix (of discretised data) with dimensions (samples x channels x epoch-trials)
			Holds data for one fly, one condition
		tau = integer - the lag between current and future states
			e.g. 1 means that the current and future sample are adjacent
			e.g. 2 means that there is one sample between the current and future samples, etc.
		n_values = number of states each *node* can be in (e.g. 2 for ON and OFF)
	Outputs:
		tpm = matrix with dimensions (n_values^channels x n_values^channels)
			Each row holds the probabilities of a past state transitioning into future states (columns)
	"""
	
	import pyphi as pyphi
	import numpy as np
	
	# Determine number of system states
	n_states = n_values ** fly_data.shape[1]
	
	# Declare TPM (all zeros)
	tpm = np.zeros((n_states, n_states))
	
	"""
	TPM Indexing (LOLI):
	e.g. for 4x4 TPM:
	
	0 = 00
	1 = 10
	2 = 01
	3 = 11
	
	Use pyphi.convert.state2loli_index(tuple) to get the index
	"""
	
	# Declare transition counter (we will divide the sum of occurrences by this to get empirical probability)
	transition_counter = np.zeros((n_states, 1))
	
	for trial in range(0, fly_data.shape[2]):
		for sample in range(0, fly_data.shape[0]-tau): # The last sample to transition is the tauth last sample (second last if tau==1) (remember that the end of range() is exclusive)
			sample_current = fly_data[sample, :, trial]
			sample_future = fly_data[sample+tau, :, trial]
			
			# Identify current state
			state_current = pyphi.convert.state2loli_index(tuple(sample_current))
			
			# Identify future state
			state_future = pyphi.convert.state2loli_index(tuple(sample_future))
			
			# Increment TPM transition by 1
			tpm[state_current, state_future] += 1
			
			# Increment transition counter
			transition_counter[state_current] += 1
	
	# Divide elements in TPM by transition counter
	# If counter is 0, then transition never occurred - to avoid dividing 0 by 0, we set the counter to 1
	counter_position = 0
	for state, counter in enumerate(transition_counter):
		if counter == 0:
			transition_counter[state] = 1
			tpm[state, :] = 1 / tpm.shape[1] # maximum entropy if no observations
		counter_position = counter_position + 1
	tpm /= transition_counter # Check vector operation
	
	return tpm# pyphi.convert.state_by_state2state_by_node(tpm)

def build_tpm_sbn(fly_data, tau, n_values):
	"""
	Builds a state-by-node TPM, holding the probabilities of each node being "on" given some past
	network states
	
	http://pyphi.readthedocs.io/en/stable/conventions.html?highlight=state%20by%20node
	
	Inputs:
		fly_data = matrix (of discretised data) with dimensions (samples x channels x epoch-trials)
			Holds data for one fly, one condition
		tau = integer - the lag between current and future states
			e.g. 1 means that the current and future sample are adjacent
			e.g. 2 means that there is one sample between the current and future samples, etc.
		n_values = number of states each *node* can be in (e.g. 2 for ON and OFF)
	Outputs:
		tpm = matrix with dimensions (n_values^channels x channels)
			Each row holds the probabilities of each channel being "on" in the future, given a past
			network state
	"""
	
	import pyphi as pyphi
	import numpy as np
	
	# Determine number of system states
	n_states = n_values ** fly_data.shape[1]
	
	# Declare TPM (all zeros)
	tpm = np.zeros((n_states, fly_data.shape[1]))
	
	"""
	TPM Indexing (LOLI):
	e.g. for 4x4 TPM:
	
	0 = 00
	1 = 10
	2 = 01
	3 = 11
	
	Use pyphi.convert.state2loli_index(tuple) to get the index (v0.8.1)
	Use pyphi.convert.state2le_index(tuple) in v1
	"""
	
	# Declare transition counter (we will divide the sum of occurrences by this to get empirical probability)
	transition_counter = np.zeros((n_states, 1))
	
	for trial in range(0, fly_data.shape[2]):
		for sample in range(0, fly_data.shape[0]-tau): # The last sample to transition is the tauth last sample (second last if tau==1) (remember that the end of range() is exclusive)
			sample_current = fly_data[sample, :, trial]
			sample_future = fly_data[sample+tau, :, trial]
			
			# Identify current state
			state_current = pyphi.convert.state2loli_index(tuple(sample_current))
			
			# Future boolean state
			sample_future_bool = sample_future.astype(bool)
			
			# Increment 'on' channels by 1
			tpm[state_current, sample_future_bool] += 1
			
			# Increment transition counter
			transition_counter[state_current] += 1
	
	# Divide elements in TPM by transition counter
	# If counter is 0, then transition never occurred - to avoid dividing 0 by 0, we set the counter to 1
	for state, counter in enumerate(transition_counter):
		if counter == 0:
			transition_counter[state] = 1
			tpm[state, :] = 1 / n_values # maximum entropy if no observations
	tpm /= transition_counter # This division works because of how we declared the vector ((n_states x 1) matrix)
	
	return tpm

def build_tpm_bin_offsets(fly_data, n_values, tau, binarise_method):
	"""
	Builds a state-by-node TPM, holding the probabilities of each node being "on" given some past
	network states. Averages across tau samples (into bins) before finding transition probabilities.
	Computes transition probabilities from each possible binning (i.e. using all offsets before binning;
	e.g. starting from sample 1, then from sample, up to tau-1)
	
	Inputs:
		fly_data = matrix (of non-discretised data) with dimensions (samples x channels x epoch-trials)
			Holds data for one fly, one condition
		tau = integer - level at which to coarse grain
			e.g. 1 means no averaging/binning across samples
			e.g. 2 means averaging each group of 2 samples
		n_values = number of states each *node* can be in (e.g. 2 for ON and OFF)
		binarise_method = string 'median' for binarising based on median or 'diff' for based on gradient
	Outputs:
		tpm = matrix with dimensions (n_values^channels x channels)
	"""
	
	import pyphi as pyphi
	import numpy as np
	
	tau_step = 1
	
	# Determine number of system states
	n_states = n_values ** fly_data.shape[1]
	
	# Declare TPM (all zeros)
	tpm = np.zeros((n_states, fly_data.shape[1]))
	
	"""
	TPM Indexing (LOLI):
	e.g. for 4x4 TPM:
	
	0 = 00
	1 = 10
	2 = 01
	3 = 11
	
	Use pyphi.convert.state2loli_index(tuple) to get the index (v0.8.1)
	Use pyphi.convert.state2le_index(tuple) in v1
	"""
	
	# Declare transition counter (we will divide the sum of occurrences by this to get empirical probability)
	transition_counter = np.zeros((n_states, 1))
	
	# For each possible sample offset before binning, bin and contribute transition probabilities
	for offset in range(0, tau):
		fly_data_offset = fly_data[offset:, :]
		
		# Downsample by averaging in bins of length tau
		fly_data_offset = tau_resample(fly_data_offset[:, :, None, None, None], tau)
		fly_data_offset = fly_data_offset[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions
		
		if binarise_method == 'median':
			# Binarise by median split
			fly_data_offset, n_values, medians = binarise_trial_median(fly_data_offset[:, :, None, None, None])
			fly_data_offset = fly_data_offset[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions
		else: # binarise_method == 'diff'
			fly_data_offset, n_values, medians = binarise_trial_diff(fly_data_offset[:, :, None, None, None])
			fly_data_offset = fly_data_offset[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions
			
		
		for sample in range(0, fly_data_offset.shape[0]-tau_step): # The last sample to transition is the second last one
			sample_current = fly_data_offset[sample, :]
			sample_future = fly_data_offset[sample+tau_step, :]
			
			# Identify current state
			state_current = pyphi.convert.state2loli_index(tuple(sample_current))
			
			# Future boolean state
			sample_future_bool = sample_future.astype(bool)
			
			# Increment 'on' channels by 1
			tpm[state_current, sample_future_bool] += 1
			
			# Increment transition counter
			transition_counter[state_current] += 1
			
	# Divide elements in TPM by transition counter
	# If counter is 0, then transition never occurred - to avoid dividing 0 by 0, we set the counter to 1
	for state, counter in enumerate(transition_counter):
		if counter == 0:
			transition_counter[state] = 1
			tpm[state, :] = 1 / n_values # maximum entropy if no observations
	tpm /= transition_counter # This division works because of how we declared the vector ((n_states x 1) matrix)
	
	return tpm, transition_counter

def build_tpm_sbn_normalise(fly_data, tau, n_values, n_normalise):
	"""
	Builds a state-by-node TPM, holding the probabilities of each node being "on" given some past
	network states
	
	http://pyphi.readthedocs.io/en/stable/conventions.html?highlight=state%20by%20node
	
	The TPM is "normalised" such that observed probabilities are corrected towards chance level as
	the number of observations supporting the probability becomes much less than n_normalise
	
	Inputs:
		fly_data = matrix (of discretised data) with dimensions (samples x channels x epoch-trials)
			Holds data for one fly, one condition
		tau = integer - the lag between current and future states
			e.g. 1 means that the current and future sample are adjacent
			e.g. 2 means that there is one sample between the current and future samples, etc.
		n_values = number of states each *node* can be in (e.g. 2 for ON and OFF)
		n_normalise = number of observations which need to be fulfilled in each transition to completely avoid
			correction/normalisation of the observed probability
	Outputs:
		tpm = matrix with dimensions (n_values^channels x channels)
			Each row holds the probabilities of each channel being "on" in the future, given a past
			network state
	"""
	
	import pyphi as pyphi
	import numpy as np
	
	# Determine number of system states
	n_states = n_values ** fly_data.shape[1]
	
	# Declare TPM (all zeros)
	tpm = np.zeros((n_states, fly_data.shape[1]))
	
	"""
	TPM Indexing (LOLI):
	e.g. for 4x4 TPM:
	
	0 = 00
	1 = 10
	2 = 01
	3 = 11
	
	Use pyphi.convert.state2loli_index(tuple) to get the index (v0.8.1)
	Use pyphi.convert.state2le_index(tuple) in v1
	"""
	
	# Declare transition counter (we will divide the sum of occurrences by this to get empirical probability)
	transition_counter = np.zeros((n_states, 1))
	
	for trial in range(0, fly_data.shape[2]):
		for sample in range(0, fly_data.shape[0]-tau): # The last sample to transition is the tauth last sample (second last if tau==1) (remember that the end of range() is exclusive)
			sample_current = fly_data[sample, :, trial]
			sample_future = fly_data[sample+tau, :, trial]
			
			# Identify current state
			state_current = pyphi.convert.state2loli_index(tuple(sample_current))
			
			# Future boolean state
			sample_future_bool = sample_future.astype(bool)
			
			# Increment 'on' channels by 1
			tpm[state_current, sample_future_bool] += 1
			
			# Increment transition counter
			transition_counter[state_current] += 1
	
	# Pad probabilities by adding transitions with probability 1/n_values (max entropy) to actual observations
	for state in range(0, tpm.shape[0]):
		pad_samples = n_normalise - transition_counter[state]
		tpm[state, :] = tpm[state, :] + (0.5*pad_samples)
		transition_counter[state] += pad_samples # After padding, number of transitions from each state should always equal n_normalise
	
	# Divide elements in TPM by transition counter
	# If counter is 0, then transition never occurred - to avoid dividing 0 by 0, we set the counter to 1
	for state, counter in enumerate(transition_counter):
		if counter == 0:
			transition_counter[state] = 1
			tpm[state, :] = 1 / n_values # Maximum entropy if no observations
	tpm /= transition_counter # This division works because of how we declared the vector ((n_states x 1) matrix)
	
	return tpm, transition_counter

def sort_constellation_by_mechanism(constellation):
	"""
	Sorts a constellation's mechanisms
	
	Mechanisms are sorted as to give an "ordered powerset order" (with the smallest sets first, etc)
	Depends on stable sorting by python
	
	Inputs:
		constellation: constellation dict from compute.pyphi.big_mip() ('partitioned_constellation' or 'unpartitioned_constellation')
	Outputs:
		sorted_const: constellation dict with sorted mechanisms
	"""
	
	# Standard sorting of lists
	sorted_const = sorted(constellation, key=lambda k: k['mechanism'])
	# Sort by length of lists (sorting is stable so previous order is maintained within second sorting)
	sorted_const = sorted(sorted_const, key=lambda k: len(k['mechanism']))
	
	return sorted_const
	

def sort_constellation_by_mechanism_orig(constellation):
	"""
	Sorts a constellation's mechanisms
	
	Mechanisms are sorted as to give an "ordered powerset order" (with the smallest sets first, etc)
	Depends on stable sorting by python
	
	Inputs:
		constellation: constellation from compute.pyphi.big_mip() (.partitioned_constellation or .unpartitioned_constellation)
	Outputs:
		sorted_const: constellation dict with sorted mechanisms
	"""
	
	# Standard sorting of lists
	sorted_const = sorted(constellation, key=lambda k: k.mechanism)
	# Sort by length of lists (sorting is stable so previous order is maintained within second sorting)
	sorted_const = sorted(sorted_const, key=lambda k: len(k.mechanism))
	
	return sorted_const

# Old functions which may no longer be useful ################################################################################

# def concat_arrays_in_array(matrix_array):
	# """
	# Concatenates each matrix (numpy array, as stored in a numpy array) into a large matrix along the 1st dimension (number of elements in the 1st dimension stays the same; concatenates the other dimensions)
	# Inputs:
		# numpy array of numpy arrays - the first dimension of these arrays should be constant
			# e.g. matrix_array[0] gives the first array
			# e.g. matrix_array[1] gives the second array, etc.
	# Outputs:
	# """
	# import numpy as np
	
	# # Take first trial
	# concatenated = matrix_array[0]
	
	# # Concatenate remaining trials
	# for matrix in range(1, matrix_array.shape[0]):
		# concatenated = np.concatenate((concatenated, matrix_array[matrix]), axis=1)
	
	# return concatenated

# # Binarising functions
# def threshold_binarise_arrays_in_array(matrix_array, thresholds):
	# """
	# Binarises matrices of an array, using the provided thresholds - each threshold corresponds to each element in the first dimension
	# Each element in the second dimension is binarised using the threshold for their position in the first dimension
	# Inputs:
		# matrix_array - array of 2D matrices, all with the same number of elements in the first dimension
		# thresholds - 1D vector which holds the thresholds to binarise by
			# Should have the same number of elements as the first dimensions of the matrices
	# Outputs:
	# """
	# import sklearn.preprocessing as skp
	# import copy
	
	# matrix_array = copy.deepcopy(matrix_array) # We don't want affect scope outside the function
	
	# for matrix in range(matrix_array.shape[0]): # i.e. for each matrix
		# for threshold_dimension in range(thresholds.shape[0]): # i.e. for each row (threshold_dimension is a D1 (row) counter)
			# matrix_array[matrix][threshold_dimension, :] = skp.binarize(matrix_array[matrix][threshold_dimension, :].reshape(1, -1), thresholds[threshold_dimension]).astype(int)
	
	# return matrix_array

# # TPM functions
# def tpm_window(matrix_array, window_start_sample, window_end_sample):
	# """
	# For a given time range (from window_start to window_end, inclusive), build a transitional probability matrix (TPM)
	# Iterates through each matrix in the array
	# Extracts values at each dimension-1 element, and at each dimension-2 element (corresponding to all channels, and the limited time window, respectively)
	# Inputs:
		# window_start_sample - which sample (column) to start building from
			# If negative, all samples are used (and window_end is ignored)
		# window_end_sample - which sample (column) to end building at (it is assumed that window_end doesn't transition to a future state)
	# Outputs:
	# """
	# import pyphi as pyphi
	# import numpy as np
	
	# # Determine number of states
	# n_states = 2**matrix_array[0].shape[0] # This assumes binarisation - a node is either ON or OFF
	
	# # Declare TPM (all zeros)
	# tpm = np.zeros((n_states, n_states))
	
	# """
	# TPM Indexing (LOLI):
	# e.g. for 4x4 TPM:
	
	# 0 = 00
	# 1 = 10
	# 2 = 01
	# 3 = 11
	
	# Use pyphi.convert.state2loli_index(tuple) to get the index
	# """
	
	# # Declare transition counter (0)
	# transition_counter = np.zeros((n_states, 1))
	
	# for matrix_counter in range(matrix_array.shape[0]): # For each matrix
		# matrix = matrix_array[matrix_counter]
		# last_sample = matrix.shape[1] - 1 # Last sample which transitions (i.e. the second last sample)
		
		# if window_start_sample == -1:
			# window_start = 0
			# window_end = last_sample - 1
		# else:
			# window_start = window_start_sample
			# window_end = window_end_sample
		
		# for sample_counter in range(window_start, window_end): # For each sample in the matrix (with a future sample)
			# sample_current = matrix[:, sample_counter] # Current state (given by the column, assuming that each row is the time series of a channel)
			# sample_future = matrix[:, sample_counter+1] # Future state
			
			# # Identify current state
			# state_current = pyphi.convert.state2loli_index(tuple(np.ndarray.tolist(sample_current)))
			
			# # Identify the following state
			# state_future = pyphi.convert.state2loli_index(tuple(np.ndarray.tolist(sample_future)))
			
			# # Increment TPM transition by 1
			# tpm[state_current, state_future] += 1
			
			# # Increment transition counter
			# transition_counter[state_current] += 1
	
	# # Divide elements in TPM by transition counter
	# tpm /= transition_counter
	
	# return pyphi.convert.state_by_state2state_by_node(tpm)

# # Phi functions
# def nchannel_phi(matrix_array, n_channels):
	# """
	# For each combination of n_channels (n_channels <= number of rows/channels in each matrix), calculates
	# the TPM across all trials using all samples, then calculates phi at each transitioning sample

	# Assumes that all channels are connected
	
	# Inputs:
		# matrix_array - array of 2D matrices, all with the same number of elements in the first dimension
		# n_channels - how many channels to calculate phi over
	# Outputs:
	# """
	# import pyphi as pyphi
	# import numpy as np
	# import itertools
	
	# combo_phis = dict()
	
	# # Get the combinations of channels
	# channels = np.arange(matrix_array[0].shape[0]) # Get list of channels
	# channel_combos = [tuple(combination) for combination in itertools.combinations(channels, n_channels)] # Get combos of channels
	
	# for combo in channel_combos:
		# print(combo)
		# # Extract relevant channels (i.e. channels specificed by the combo) across all trials
		# samples_all = matrix_array[0][combo, :]
		# for matrix_counter in range(1, matrix_array.size):
			# np.concatenate((samples_all, matrix_array[matrix_counter][combo, :]), axis=1)
		# # Build the TPM from those channels
		# tpm_input = np.empty((1), dtype=object)
		# tpm_input[0] = samples_all
		# tpm = tpm_window(tpm_input, -1, 0) # We place the concatenated matrix into an array as input for tpm_window(), which takes an array of matrices
		
		# #print("computing phi")
		# # Compute phi at each trial sample
		# network = pyphi.Network(tpm)
		# phis = np.empty((matrix_array.size), dtype=object)
		# for trial_counter in range(1, 2):#matrix_array.size):
			# phis[trial_counter] = np.empty(matrix_array[trial_counter].shape[1])
			# trial = matrix_array[trial_counter]
			# for sample_counter in range(matrix_array[trial].shape[1]):
				# print("trial" + str(trial_counter) + " sample" + str(sample_counter))
				# sample = trial[combo, sample_counter]
				# #print(sample)
				# state = tuple(np.ndarray.tolist(sample))
				# #print("state done")
				# subsystem = pyphi.Subsystem(network, state, range(network.size))
				# #print("subsystem done")
				# phis[trial_counter][sample_counter] = pyphi.compute.big_phi(subsystem)
				# #print("phi = " + str(phis[trial_counter][sample_counter]))
		
		# # Add combo results to dictionary
		# combo_phis[combo] = phis
	
	# return combo_phis

# # Plotting functions
# def plot_mean(matrix):
	# """
	# Calculates the mean and standard deviation across the first dimension (axis 0), and plots across the second dimension (axis 1)
	# Inputs: numpy 2D array
	# """
	# import numpy as np
	# import matplotlib.pyplot as plt
	
	# if len(matrix.shape) > 2:
		# return
	
	# if len(matrix.shape) == 2: # two dimensions provided
		# x_values = np.arange(0, matrix.shape[1])
		# mean_values = matrix.mean(axis=0)
		# std_values = matrix.std(axis=0)
	# else: # only one dimension was provided
		# x_values = np.arange(0, matrix.shape[0])
		# mean_values = matrix
		# std_values = matrix
	
	# plt.plot(x_values, mean_values)
	# if len(matrix.shape) == 2: # two dimensions provided
		# plt.fill_between(x_values, mean_values - std_values, mean_values + std_values, alpha=0.5, linewidth=0)
	# plt.autoscale(enable=True, axis='x', tight=True)
	# plt.show(block=False)


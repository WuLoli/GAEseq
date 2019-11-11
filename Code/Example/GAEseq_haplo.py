import tensorflow as tf
import numpy as np
from Bio import SeqIO
import random
from itertools import permutations
from helper import *
import time
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

# import config file
with open('config', 'r') as f:
	df = f.readlines()

data = []
for item in df:
	data.append(item.strip('\n').split(' : '))

zone_name = data[8][1] # zone name
SNVmatrix_name = zone_name + '_SNV_matrix.txt' # SNP matrix file name

# import SNP fragment matrix
with open(SNVmatrix_name, 'r') as f:
	SNVmatrix_list = f.readlines()

SNVmatrix = np.zeros((len(SNVmatrix_list), len(SNVmatrix_list[0][0:-1:2])))

for i in range(len(SNVmatrix_list)):
	SNVmatrix[i, :] = np.fromstring(SNVmatrix_list[i][0:-1:1], dtype = int, sep = ' ')

# observed matrix
observed_matrix = SNVmatrix.copy()
projection_matrix = (observed_matrix != 0) * 1

# transfer SNPmatrix to unfolded tensor
observed_tensor = OneHotEncode(SNVmatrix)
projection_tensor = np.tile(projection_matrix, (1, 4))

# number of different nucleotides
num_allele = 4

# training parameters
learning_rate = 0.0001 # step size
training_epoch = 100 # number of epochs
num_read, len_haplo = observed_matrix.shape # number of reads; number of nucleotides
tiny_value = 10 ** -12 # a very small value

# adjacency matrix
adjacency_matrix = np.zeros((observed_matrix.shape[0], observed_matrix.shape[1], num_allele))
for i in range(num_allele):
	adjacency_matrix[:, :, i] = (observed_matrix == i + 1) * 1

# factor
read2haplo_factor = np.zeros((len_haplo, num_read, num_allele))
for i in range(num_allele):
	for j in range(len_haplo):
		read2haplo_factor[j, :, i] = (adjacency_matrix[:, :, i].T)[j, :] / np.sum(observed_matrix[:, j] != 0)

haplo2read_factor = np.zeros((num_read, len_haplo, num_allele))
for i in range(num_allele):
	for j in range(num_read):
		haplo2read_factor[j, :, i] = adjacency_matrix[j, :, i] / np.sum(observed_matrix[j, :] != 0)

ploidy = int(data[9][1]) # ploidy


# dimension of each layer
num_layers = int(data[10][1])
C_i = [int(len_haplo - i / (num_layers + 1) * (len_haplo - ploidy)) for i in range(1, num_layers + 1)]

# Tensorflow Graph input
tf.reset_default_graph()

# Tensorflow placeholder
X = tf.placeholder(tf.float32, [num_read, len_haplo], name = 'X') # placeholder for SNP fragment matrix
V = tf.placeholder(tf.float32, [ploidy, len_haplo * 4], name = 'V') # placeholder for haplotype matrix in unfoldered tensor structure
keep_prob_u2i = tf.placeholder(tf.float32, name = 'keep_prob_u2i') # dropout probability of read nodes layer to SNP nodes layer
keep_prob_i2u = tf.placeholder(tf.float32, name = 'keep_prob_i2u') # dropout probability of SNP nodes layer to read nodes layer
keep_prob_dense = tf.placeholder(tf.float32, name = 'keep_prob_dense') # dropout probability of dense layer

# define weights and bias for each layer 
Wr_item_1 = tf.get_variable('Wr_item_1', shape = [len_haplo, C_i[0], num_allele], initializer = tf.contrib.layers.xavier_initializer())
Br_item_1 = tf.get_variable('Br_item_1', shape = [len_haplo, C_i[0], num_allele], initializer = tf.contrib.layers.xavier_initializer())

for i in range(2, num_layers + 1):
	variable_name = 'Wr_item_' + str(i)
	exec(variable_name + ' = tf.get_variable(' + '\'' + variable_name + '\''  + ', shape = [C_i[' + str(i - 2) + '], C_i[' + str(i - 1) + '], num_allele], initializer = tf.contrib.layers.xavier_initializer())')	
	variable_name = 'Br_item_' + str(i)
	if i % 2 == 0:
		exec(variable_name + ' = tf.get_variable(' + '\'' + variable_name + '\''  + ', shape = [num_read, C_i[' + str(i - 1) + '], num_allele], initializer = tf.contrib.layers.xavier_initializer())')		 
	else:
		exec(variable_name + ' = tf.get_variable(' + '\'' + variable_name + '\''  + ', shape = [len_haplo, C_i[' + str(i - 1) + '], num_allele], initializer = tf.contrib.layers.xavier_initializer())')		 
	
W_dense = tf.get_variable('W_dense', shape = [C_i[-1], ploidy], initializer = tf.contrib.layers.xavier_initializer())

# define read nodes to SNP nodes layer
def read2haplo_layer(W, B, x, C_before, C_after, num_allele = num_allele, read2haplo_factor = read2haplo_factor):
	res = tf.Variable(tf.zeros(shape = [len_haplo, C_after]))
	
	for i in range(num_allele):
		tem_v = tf.matmul(tf.cast(read2haplo_factor[:, :, i], tf.float32), x)
		tem_v = tf.matmul(tem_v, W[:, :, i])
		tem_v = tf.math.add(tem_v, B[:, :, i])
		tem_v = tf.nn.dropout(tf.nn.relu(tem_v), keep_prob_u2i)
		res += tem_v
	
	return res

# define SNP nodes to read nodes layer
def haplo2read_layer(W, B, x, C_before, C_after, num_allele = num_allele, haplo2read_factor = haplo2read_factor):
	res = tf.Variable(tf.zeros(shape = [num_read, C_after]))
	
	for i in range(num_allele):
		tem_u = tf.matmul(tf.cast(haplo2read_factor[:, :, i], tf.float32), x)
		tem_u = tf.matmul(tem_u, W[:, :, i])
		tem_u = tf.math.add(tem_u, B[:, :, i])
		tem_u = tf.nn.dropout(tf.nn.relu(tem_u), keep_prob_i2u)
		res += tem_u
	
	return res	

# define graph convolutional encoder
def graph_encoder(observed_matrix, num_allele = num_allele, num_read = num_read, len_haplo = len_haplo):		
	# item layer
	M_i_1 = read2haplo_layer(Wr_item_1, Br_item_1, observed_matrix, len_haplo, C_i[0])
	
	if num_layers > 2:
		for i in range(2, num_layers):
			M_name_previous = 'M_i_' + str(i - 1)
			M_name = 'M_i_' + str(i)
			variable_name = 'Wr_item_' + str(i)
			bias_name = 'Br_item_' + str(i)
			
			if i % 2 == 0:
				exec(M_name + ' = haplo2read_layer(' + variable_name + ', ' + bias_name + ', ' + M_name_previous + ', C_i[' + str(i - 2) + '], C_i[' + str(i - 1) + '])')
			else:
				exec(M_name + ' = read2haplo_layer(' + variable_name + ', ' + bias_name + ', ' + M_name_previous + ', C_i[' + str(i - 2) + '], C_i[' + str(i - 1) + '])')
	else:
		i = 1
	
	i += 1
	M_name_previous = 'M_i_' + str(i - 1)
	variable_name = 'Wr_item_' + str(i)
	bias_name = 'Br_item_' + str(i)
	
	Hv = eval('haplo2read_layer(' + variable_name + ', ' + bias_name + ', ' + M_name_previous + ', C_i[' + str(i - 2) + '], C_i[' + str(i - 1) + '])')

	U = tf.exp(200 * tf.nn.relu(tf.nn.dropout(tf.matmul(tf.nn.relu(Hv), W_dense), keep_prob_dense)))

	sum_U = tf.stack([tf.reduce_sum(U, 1)] * W_dense.shape[1], axis = 0)
	
	return tf.divide(U, tf.transpose(sum_U))

# define haplotype decoder
def U2V(U, M_E, R = ploidy, hap_len = len_haplo):
	min_index = np.argmax(U, axis = 1)
	V_major = np.zeros((R, hap_len)) # majority voting result
	ACGTcount = ACGT_count(M_E)
	
	for i in range(R):		 
		reads_single = M_E[min_index == i, :] # all reads from one haplotypes
		single_sta = np.zeros((hap_len, 4))
		
		if len(reads_single) != 0:
			single_sta = ACGT_count(reads_single) # ACGT statistics of a single nucleotide position
		V_major[i, :] = np.argmax(single_sta, axis = 1) + 1.0			  

		uncov_pos = np.where(np.sum(single_sta, axis = 1) == 0)[0]

		for j in range(len(uncov_pos)):
			if len(np.where(ACGTcount[uncov_pos[j], :] == max(ACGTcount[uncov_pos[j], :]))[0]) != 1: # if not covered, select the most doninant one based on 'ACGTcount'	 
				tem = np.where(ACGTcount[uncov_pos[j], :] == max(ACGTcount[uncov_pos[j], :]))[0]
				V_major[i, uncov_pos[j]] = tem[int(np.floor(random.random() * len(tem)))] + 1
			else:
				V_major[i, uncov_pos[j]] = np.argmax(ACGTcount[uncov_pos[j], :]) + 1

	return V_major

# approximated read origin indicator matrix
U = graph_encoder(X, num_allele = num_allele, num_read = num_read, len_haplo = len_haplo)
U_test = tf.identity(U, name = 'U')

# define loss and optimizer
loss = 0.5 * tf.reduce_sum(tf.multiply(tf.square(tf.subtract(tf.matmul(U, V), observed_tensor)), tf.cast(projection_tensor, tf.float32)))
optimizer = tf.train.AdamOptimizer(learning_rate = learning_rate).minimize(loss)

# initialize the variables
init = tf.global_variables_initializer()

# number of experiments
num_exp = int(data[11][1])

# whether to use GPU
GPU_flag = int(data[15][1])

# keep probability
p_r2s = 1 - int(data[12][1])
p_r2s = 1 - int(data[13][1])
p_dense = 1 - int(data[14][1])

# start training auto-encoder
MEC_experiment = []
haplo_experiment = np.zeros((ploidy, len_haplo, num_exp))
flag = 10 ** 8
time_record = []

for iteration in range(num_exp):
	start_time = time.time()
	saver = tf.train.Saver(max_to_keep = training_epoch)

	if GPU_flag == -1:
		setting = tf.ConfigProto(allow_soft_placement = True, log_device_placement = False)
	else:
		gpu_options = tf.GPUOptions(visible_device_list = str(GPU_flag))
		setting = tf.ConfigProto(allow_soft_placement = True, log_device_placement = False, gpu_options = gpu_options)

	with tf.Session(config = setting) as sess:
		
		sess.run(init)
		
		index_matrix = sess.run(U, feed_dict = {X: observed_matrix, keep_prob_u2i: p_r2s, keep_prob_i2u: p_r2s, keep_prob_dense : p_dense})
		haplo_matrix = U2V(index_matrix, observed_matrix) 

		epoch_haplo = np.zeros((ploidy, len_haplo, training_epoch))
		epoch_MEC = []

		for epoch in range(training_epoch):
			sess.run(optimizer, feed_dict = {X: observed_matrix, V: OneHotEncode(haplo_matrix), keep_prob_u2i: 1, keep_prob_i2u: 1, keep_prob_dense : 1})
			index_matrix = sess.run(U, feed_dict = {X: observed_matrix, keep_prob_u2i: 1, keep_prob_i2u: 1, keep_prob_dense : 1})
			haplo_matrix = U2V(index_matrix, observed_matrix)  
			epoch_haplo[:, :, epoch] = haplo_matrix
			epoch_MEC.append(MEC(SNVmatrix, haplo_matrix))

	# calibration
	try:
		GAE_hap = epoch_haplo[:, :, np.argmin(np.array(epoch_MEC))]
		index = []
		for i in range(SNVmatrix.shape[0]):
			dis = np.zeros((GAE_hap.shape[0]))
			for j in range(GAE_hap.shape[0]):
				dis[j] = hamming_distance(SNVmatrix[i, :], GAE_hap[j, :])
			index.append(np.argmin(dis))

		new_haplo = np.zeros((GAE_hap.shape))
		for i in range(GAE_hap.shape[0]):
			new_haplo[i, :] = np.argmax(ACGT_count(SNVmatrix[np.array(index) == i, :]), axis = 1) + 1

		MEC_experiment.append(MEC(SNVmatrix, new_haplo))
		haplo_experiment[:, :, iteration] = new_haplo
	except:
		MEC_experiment.append(MEC(SNVmatrix, epoch_haplo[:, :, np.argmin(np.array(epoch_MEC))]))
		haplo_experiment[:, :, iteration] = epoch_haplo[:, :, np.argmin(np.array(epoch_MEC))]

	# # record MEC results
	# with open(zone_name + '_MEC_result_{}_Strains.txt'.format(ploidy), 'a') as f:
	# 	f.write('Experiment' + str(iteration + 1) + ' : ' + str(MEC_experiment[-1])  + '\n')

	end_time = time.time()
	time_record.append(end_time - start_time)
	print('Checking Rank: {} - '.format(ploidy) + 'Experiment: {}/{} - ETA: {}s '.format(iteration + 1, num_exp, int(np.mean(time_record) * (num_exp - (iteration + 1)))), end = '\r')

result = haplo_experiment[:, :, np.argmin(np.array(MEC_experiment))]

with open(zone_name + '_Haplo_{}_Strains.txt'.format(ploidy), 'w') as f:
	for i in range(result.shape[0]):
		f.write('Haplotype ' + str(i + 1) + '\n')
		for j in range(result.shape[1]):
			if result[i, j] == 1:
				f.write('A')
			elif result[i, j] == 2:
				f.write('C')
			elif result[i, j] == 3:
				f.write('G')
			elif result[i, j] == 4:
				f.write('T')
			elif result[i, j] == 0:
				f.write('*')
		f.write('\n')
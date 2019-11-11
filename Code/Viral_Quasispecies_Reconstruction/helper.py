import tensorflow as tf
import numpy as np
from Bio import SeqIO
import random
from itertools import permutations

# convert list of sequence to numpy array (A -> 1, C -> 2, G -> 3, T -> 4 and unobserved -> 0)
def list2array(ViralSeq_list):
    ViralSeq = np.zeros((len(ViralSeq_list), len(ViralSeq_list[0])))

    for i in range(len(ViralSeq_list)):
        for j in range(len(ViralSeq_list[0])):
            if ViralSeq_list[i][j] == 'A':
                ViralSeq[i, j] = 1
            elif ViralSeq_list[i][j] == 'C':
                ViralSeq[i, j] = 2
            elif ViralSeq_list[i][j] == 'G':
                ViralSeq[i, j] = 3
            elif ViralSeq_list[i][j] == 'T':
                ViralSeq[i, j] = 4
    
    return ViralSeq

# calculate hamming distance between two sequences
def hamming_distance(read, haplo):
    return sum((haplo - read)[np.where(read != 0)] != 0)

# evaluate the recall rate and reconstruction rate
def recall_reconstruction_rate(Recovered_Haplo, SNVHaplo):
    distance_table = np.zeros((len(Recovered_Haplo), len(SNVHaplo)))
    for i in range(len(Recovered_Haplo)):
        for j in range(len(SNVHaplo)):
            distance_table[i, j] = hamming_distance(SNVHaplo[j, :], Recovered_Haplo[i, :])
    
    index = list(permutations(list(range(SNVHaplo.shape[0]))))
    distance = []
    for item in index:
        count = 0
        for i in range(len(item)):
            count += distance_table[i, item[i]]
        distance.append(count)
    index = index[np.argmin(np.array(distance))]
    
    reconstruction_rate = []
    for i in range(len(index)):
        reconstruction_rate.append(1 - distance_table[i, index[i]] / SNVHaplo.shape[1])
    
    recall_rate = reconstruction_rate.count(1) / len(reconstruction_rate)
    
    CPR = 1 - min(distance) / (len(distance_table) * SNVHaplo.shape[1])
    
    return distance_table, reconstruction_rate, recall_rate, index, CPR

# one hot encode
def OneHotEncode(Haplo):
    return np.concatenate((np.double(Haplo == 1), np.double(Haplo == 2), np.double(Haplo == 3), np.double(Haplo == 4)), axis = 1)

# one hot decode
def OneHotDecode(Haplo):
    res = np.zeros((len(Haplo), int(len(Haplo[0]) / 4)))
    Haplo = Haplo.reshape(len(Haplo), int(len(Haplo[0]) / 4), 4, order = 'F')
    for i in range(Haplo.shape[0]):
        for j in range(Haplo.shape[1]):
            if max(Haplo[i, j, :]) != 0:
                res[i, j] = np.argmax(Haplo[i, j, :]) + 1      
    return res

# get the ACGT statistics of a read matrix
def ACGT_count(M_E):
    out = np.zeros((len(M_E[0, :]), 4))
    for i in range(4):
        out[:, i] = (M_E == (i + 1)).sum(axis = 0)

    return out 

# use haplotype matrix to recover observed matrix
def haplo2complete(observed_matrix, haplo_matrix):
    res = np.zeros((observed_matrix.shape))
    
    for i in range(observed_matrix.shape[0]):
        distance = np.zeros((1, haplo_matrix.shape[0]))
        for j in range(haplo_matrix.shape[0]):
            distance[0, j] = hamming_distance(observed_matrix[i, :], haplo_matrix[j, :])
        res[i, :] = haplo_matrix[np.argmin(distance), :]
        
    return res        

# calculate MEC
def MEC(SNVmatrix, Recovered_Haplo):
    res = 0
    
    for i in range(len(SNVmatrix)):
        dis = [hamming_distance(SNVmatrix[i, :], Recovered_Haplo[j, :]) for j in range(len(Recovered_Haplo))]
        res += min(dis)
        
    return res
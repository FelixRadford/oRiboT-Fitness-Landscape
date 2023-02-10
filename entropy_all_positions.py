import pandas as pd

import numpy as np


df_matrix = pd.read_csv ('Library_Sites_Matrix.csv', header=None)   #reads CSV to obtain the sites matrix for this library


df_matrix[df_matrix<0] = 0      #set all negative values to 0 in this matrix
 
df2 = df_matrix.sum()            # create matrix with sums of all columns

prob_matrix = df_matrix/df2       #probabilities of all positions



prob_1 = prob_matrix[0]           #probability of 1st base

prob_1_n = prob_1[prob_1>0]      # remove zeroes

motif_1 = prob_1_n.values.tolist()    #create list of probability values for calculation of entropy

 
prob_2 = prob_matrix[1]           #probability of 1st base

prob_2_n = prob_2[prob_2>0]      # remove zeroes

motif_2 = prob_2_n.values.tolist()    #create list of probability values for calculation of entropy



prob_3 = prob_matrix[2]           #probability of 1st base

prob_3_n = prob_3[prob_3>0]      # remove zeroes

motif_3 = prob_3_n.values.tolist()    #create list of probability values for calculation of entropy



prob_4 = prob_matrix[3]           #probability of 1st base

prob_4_n = prob_4[prob_4>0]      # remove zeroes

motif_4 = prob_4_n.values.tolist()    #create list of probability values for calculation of entropy



prob_5 = prob_matrix[4]           #probability of 1st base

prob_5_n = prob_5[prob_5>0]      # remove zeroes

motif_5 = prob_5_n.values.tolist()    #create list of probability values for calculation of entropy



entropies = []         #list of entropies computed for this site, in order from position 1 to 5

def compute_entropy(motif):
    arr = np.array(motif)
    H = (arr.flatten() * np.log2(arr.flatten())).sum(axis=0) # Formula from: http://en.wikipedia.org/wiki/Sequence_logo
    entropies.append(-H)
    print 'entropy:', -H.mean(), 'bits'

#comput all entropies in bits

p_1 = compute_entropy(motif_1)
p_2 = compute_entropy(motif_2)
p_3 = compute_entropy(motif_3)
p_4 = compute_entropy(motif_4)
p_5 = compute_entropy(motif_5)





# output list to .csv file
entropies = pd.DataFrame(entropies)        #convert list to dataframe
entropies.to_csv('entropies.csv', index=False, header=None)        #write to file



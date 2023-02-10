
import csv
import pandas as pd
import os
import subprocess

#Analyze positive, negative, and neutral portions of the library

process1 = subprocess.Popen(["python", "CorMatrix_large.py"], cwd='Positive/') # Create and launch process pop.py using python interpreter; specify current working directory
process2 = subprocess.Popen(["python", "CorMatrix_large.py"], cwd='Negative/')
process3 = subprocess.Popen(["python", "CorMatrix_large.py"], cwd='Neutral/')

process1.wait() # Wait for process1 to finish (basically wait for script to finish)
process2.wait()
process3.wait()

# put positives into pd.dataframe from csv file in "Positives/"

df_pos = pd.read_csv ('Positive/cormatrix.csv', header=None) #reads CSV for positives matrix, first row starts at 0

Pos_Sum = df_pos.values.sum() #sum of all values in positive matrix

Pos_norm1 = df_pos.div(Pos_Sum)   #divide positives by sum to create new table



# Find sum of negatives and divide negatives by sum to create new table

df_neg = pd.read_csv ('Negative/cormatrix.csv', header=None) #reads CSV for positives matrix, first row starts at 0

Neg_Sum = df_neg.values.sum() #sum of all values in positive matrix

Neg_norm1 = df_neg.div(Neg_Sum)   #divide positives by sum to create new table

 


# Find sum of neutrals and divide neutrals by sum to create new table

df_n = pd.read_csv ('Neutral/cormatrix.csv', header=None) #reads CSV for positives matrix, first row starts at 0

n_Sum = df_n.values.sum() #sum of all values in positive matrix

n_norm1 = df_n.div(n_Sum)   #divide positives by sum to create new table

 

# Normalize positives by neutrals

Pos_norm2 = Pos_norm1.div(n_norm1)
 



# Normalize negatives by neutrals

Neg_norm2 = Neg_norm1.div(n_norm1)

# Divide positives by negatives to generate final matrix

Normlz = Pos_norm2.sub(Neg_norm2)    #subtract normalized negatives matrix from normalized positives matrix
Comb_data = Normlz.mul(100)          #multiply subtracted by 100 to obtain percantage


# Write final matrix to CSV

Comb_data.to_csv('Library_covariation_matrix.csv', index=False, header=None)


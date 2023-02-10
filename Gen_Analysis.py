
import csv
import pandas as pd
import os
import subprocess

# Define the directory where file to analyze is located

directory = '' #specify directory where .fastq file for analysis is located

#The address of the directory where the file to be analyzed is located is written to "address.txt", which is then read by each of the three analysis scripts
directory_file1 = open("Positive/address.txt", "w")
directory_file1.write(directory)
directory_file1.close()

directory_file2 = open("Negative/address.txt", "w")
directory_file2.write(directory)
directory_file2.close()

directory_file3 = open("Neutral/address.txt", "w")
directory_file3.write(directory)
directory_file3.close()


#Analyze positive, negative, and neutral portions of the library

process1 = subprocess.Popen(["python", "Positive_Analysis.py"], cwd='Positive/') # Create and launch process pop.py using python interpreter; specify current working directory
process2 = subprocess.Popen(["python", "Negative_Analysis.py"], cwd='Negative/')
process3 = subprocess.Popen(["python", "Neutral_Analysis.py"], cwd='Neutral/')

process1.wait() # Wait for process1 to finish (basically wait for script to finish)
process2.wait()
process3.wait()

# put positives into pd.dataframe from csv file in "Positives/"

df_pos = pd.read_csv ('Positive/Sites-Matrix.csv', header=None) #reads CSV for positives matrix, first row starts at 0

Pos_Sum = df_pos.values.sum() #sum of all values in positive matrix

Pos_norm1 = df_pos.div(Pos_Sum)   #divide positives by sum to create new table



# Find sum of negatives and divide negatives by sum to create new table

df_neg = pd.read_csv ('Negative/Sites-Matrix.csv', header=None) #reads CSV for positives matrix, first row starts at 0

Neg_Sum = df_neg.values.sum() #sum of all values in positive matrix

Neg_norm1 = df_neg.div(Neg_Sum)   #divide positives by sum to create new table

 


# Find sum of neutrals and divide neutrals by sum to create new table

df_n = pd.read_csv ('Neutral/Sites-Matrix.csv', header=None) #reads CSV for positives matrix, first row starts at 0

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

Comb_data.to_csv('Library_Sites_Matrix.csv', index=False, header=None)


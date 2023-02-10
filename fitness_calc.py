'''Calculation of fitness for all unique sequences in a library'''

import itertools
import pandas as pd
import numpy as np
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from scipy.stats import tstd
from scipy.stats import trim_mean
 


# open csv file with all sequences and their frequencies. Put positives into pd.dataframe from csv file in "Positives/"

df_pos = pd.read_csv ('Positive/dict.csv', header=None) #reads CSV for positives matrix, first row starts at 0

#neutral

df_neutr = pd.read_csv ('Neutral/dict.csv', header=None)

#negative

df_neg = pd.read_csv ('Negative/dict.csv', header=None)


# Normalize positives by neutrals
Pos_norm = df_pos
Pos_norm[1] = Pos_norm[1].sub(df_neutr[1])




# Normalize negatives by neutrals
Neg_norm = df_neg
Neg_norm[1] = Neg_norm[1].sub(df_neutr[1])



# Obtain fitness by subtracting negatives from neutral and multiplying by 100

Mutants_fitness_matrix = Pos_norm
Mutants_fitness_matrix[1]  = Pos_norm[1].div(Neg_norm[1])    #subtract negative matrix from positive matrix


#drop all infinite values
Mutants_fitness_matrix.replace([np.inf, -np.inf], np.nan, inplace=True) #replace all infinite values with nan
Mutants_fitness_matrix.dropna(inplace=True)   #drop all inf values


#drop all values from dataframe that lower than -100 and greater than 100

Mutants_fitness_matrix = Mutants_fitness_matrix.loc[(Mutants_fitness_matrix[1] > - 100)]
Mutants_fitness_matrix = Mutants_fitness_matrix.loc[(Mutants_fitness_matrix[1] < 100)]


# get mean and SD of dataset

mean = scipy.mean(Mutants_fitness_matrix[1])
SD = tstd(Mutants_fitness_matrix[1])
t_mean = trim_mean(Mutants_fitness_matrix[1], 0.05) # Trim 5% at both ends; not used for normalization

# maximum value:
my_max = Mutants_fitness_matrix[1].loc[Mutants_fitness_matrix[1].idxmax()]   


#Normalize data
Norm_data =  Mutants_fitness_matrix
Norm_data[1] =  Mutants_fitness_matrix[1].sub(mean)       
Norm_data[1] =  Mutants_fitness_matrix[1].div(SD)


#plot smoothed histogram of the data
ax = sns.kdeplot(Norm_data[1], shade=True, color="b")
#ax = sns.distplot(Norm_data[1])

plt.legend([],[], frameon=False)
plt.xlim(-14, 14)     #specificy the range shown on the graph for the x axis
plt.savefig('fitness_hist_norm.png', dpi=300)
 








##This script normalizes the covariation matrices for each library

import pandas as pd
import numpy as np
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from scipy.stats import tstd
from scipy.stats import trim_mean



# put positives into pd.dataframe from csv file in "Positives/"

df_pos = pd.read_csv ('Positive/cormatrix.csv', header=None) #reads CSV for positives matrix, first row starts at 0




# Find sum of negatives and divide negatives by sum to create new table

df_neg = pd.read_csv ('Negative/cormatrix.csv', header=None) #reads CSV for positives matrix, first row starts at 0




# Find sum of neutrals and divide neutrals by sum to create new table

df_n = pd.read_csv ('Neutral/cormatrix.csv', header=None) #reads CSV for positives matrix, first row starts at 0



# Normalize positives by neutrals

Pos_norm2 = df_pos.sub(df_n)
 



# Normalize negatives by neutrals

Neg_norm2 = df_neg.sub(df_n)

# Divide positives by negatives to generate final matrix

Normlz = Pos_norm2.div(Neg_norm2)    #subtract normalized negatives matrix from normalized positives matrix

#obtain means and SDs of each column

mean1 = scipy.mean(Normlz[0])
mean2 = scipy.mean(Normlz[1])
mean3 = scipy.mean(Normlz[2])
mean4 = scipy.mean(Normlz[3])
mean5 = scipy.mean(Normlz[4])
mean6 = scipy.mean(Normlz[5])
mean7 = scipy.mean(Normlz[6])
mean8 = scipy.mean(Normlz[7])
mean9 = scipy.mean(Normlz[8])
mean10 = scipy.mean(Normlz[9])

SD1 = tstd(Normlz[0])
SD2 = tstd(Normlz[1])
SD3 = tstd(Normlz[2])
SD4 = tstd(Normlz[3])
SD5 = tstd(Normlz[4])
SD6 = tstd(Normlz[5])
SD7 = tstd(Normlz[6])
SD8 = tstd(Normlz[7])
SD9 = tstd(Normlz[8])
SD10 = tstd(Normlz[9])



#Normalize data
Normlz[0] =  Normlz[0].sub(mean1)       
Normlz[0] =  Normlz[0].div(SD1)

Normlz[1] =  Normlz[1].sub(mean2)       
Normlz[1] =  Normlz[1].div(SD2)

Normlz[2] =  Normlz[2].sub(mean3)       
Normlz[2] =  Normlz[2].div(SD3)

Normlz[3] =  Normlz[3].sub(mean4)       
Normlz[3] =  Normlz[3].div(SD4)

Normlz[4] =  Normlz[4].sub(mean5)       
Normlz[4] =  Normlz[4].div(SD5)

Normlz[5] =  Normlz[5].sub(mean6)       
Normlz[5] =  Normlz[5].div(SD5)

Normlz[6] =  Normlz[6].sub(mean7)       
Normlz[6] =  Normlz[6].div(SD5)

Normlz[7] =  Normlz[7].sub(mean8)       
Normlz[7] =  Normlz[7].div(SD5)

Normlz[8] =  Normlz[8].sub(mean9)       
Normlz[8] =  Normlz[8].div(SD5)

Normlz[9] =  Normlz[9].sub(mean10)       
Normlz[9] =  Normlz[9].div(SD5)

# Write final matrix to CSV

Normlz.to_csv('cormatrix_normalized.csv', index=False, header=None)


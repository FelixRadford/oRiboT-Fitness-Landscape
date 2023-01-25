# Generate a matrix of pairwise epistasis for all fAB, fA, and fB values
# Felix Radford, 2022

import itertools
import pandas as pd
import csv


WT = "AAAGA"              # input WT sequence for this PTC region




# open csv file with all sequences and their frequencies, from positive, neutral, and negative selections.

#positive 

df_pos = pd.read_csv ('Positive/dict.csv', header=None) #reads CSV for positives matrix, first row starts at 0

#neutral

df_neutr = pd.read_csv ('Neutral/dict.csv', header=None)

#negative

df_neg = pd.read_csv ('Negative/dict.csv', header=None)


# Normalize positives by neutrals
Pos_norm = df_pos
Pos_norm[1] = Pos_norm[1].sub(df_neutr[1])

#Pos_norm.to_csv('Pos_norm_matrix.csv', index=False, header=None)


# Normalize negatives by neutrals
Neg_norm = df_neg
Neg_norm[1] = Neg_norm[1].sub(df_neutr[1])



# Obtain fitness by subtracting negatives from neutral and multiplying by 100

Mutants_fitness_matrix = Pos_norm
Mutants_fitness_matrix[1]  = Pos_norm[1].div(Neg_norm[1])    #subtract negative matrix from positive matrix




Mutants_fitness_matrix.to_csv('A_Mutants_fitness_matrix_all.csv', index=False, header=None)







# Lists to fill with fitness values for epistasis for each combinations of sites:

e1_2 = []

e1_3 = []

e1_4 = []

e1_5 = []

e2_3 = []

e2_4 = []

e2_5 = []

e3_4 = []

e3_5 = []

e4_5 = []

s = "n"
 



########## Make a list of all combinations of dual site mutants


nucleotide_list = ['A', 'T', 'G', 'C']

randomizer = itertools.product(nucleotide_list, repeat = 2)    #produce complete set of random combinations of two nucleotides
randomizer_list = [''.join(i) for i in randomizer]

def Ep_e1_2():
	for i in randomizer_list:
		AB = i + WT[2:5]         
		if AB == WT:
			e1_2.append(s)
		elif AB[0] == WT[0]:
			e1_2.append(s)
		elif AB[1] == WT[1]:
			e1_2.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe
			f_AB = float(dual_fitness[1].values[0])       # value of fitness corresponding to the sequence found and outputted as float
			A = i[:1] + WT[1:5] #also get fitness
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B = WT[0] + i[1:] + WT[2:5] 
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0:     #does the sequence exist?
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0])
				epistasis = f_AB - (f_A*f_B)
				e1_2.append(epistasis)

					 

def Ep_e1_3():
	for i in randomizer_list:
		AB = i[:1] + WT[1] + i[1:] + WT[3:5]      
		if AB == WT:
			e1_3.append(s)
		elif AB[0] == WT[0]:
			e1_3.append(s)
		elif AB[2] == WT[2]:
			e1_3.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe
			f_AB = float(dual_fitness[1].values[0])       # value of fitness corresponding to the sequence found and outputted as float
			A = i[:1] + WT[1:5]
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B = WT[0:2] + i[1:] +  WT[3:5]
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0:
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0])
				epistasis = f_AB - (f_A*f_B)
				e1_3.append(epistasis)

def Ep_e1_4():
	for i in randomizer_list:
		AB = i[:1] + WT[1:3] + i[1:] + WT[4]      
		if AB == WT:
			e1_4.append(s)
		elif AB[0] == WT[0]:
			e1_4.append(s)
		elif AB[3] == WT[3]:
			e1_4.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe
			f_AB = float(dual_fitness[1].values[0])       # value of fitness corresponding to the sequence found and outputted as float
			A = i[:1] + WT[1:5]
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B = WT[0:3] + i[1:] +  WT[4]
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0: 
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0])
				epistasis = f_AB - (f_A*f_B)
				e1_4.append(epistasis)

def Ep_e1_5():
	for i in randomizer_list:
		AB = i[:1] + WT[1:4] + i[1:]     
		if AB == WT:
			e1_5.append(s)
		elif AB[0] == WT[0]:
			e1_5.append(s)
		elif AB[4] == WT[4]:
			e1_5.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe     
			A = i[:1] + WT[1:5]
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B = WT[0:4] + i[1:]
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0: 
				f_AB = float(dual_fitness[1].values[0])  
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0]) 
				epistasis = f_AB - (f_A*f_B)
				e1_5.append(epistasis)


def Ep_e2_3():
	for i in randomizer_list:
		AB = WT[0] + i + WT[3:5]      
		if AB == WT:
			e2_3.append(s)
		elif AB[1] == WT[1]:
			e2_3.append(s)
		elif AB[2] == WT[2]:
			e2_3.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe
			f_AB = float(dual_fitness[1].values[0])       # value of fitness corresponding to the sequence found and outputted as float
			A = WT[0] + i[:1] + WT[2:5]
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B = WT[0:2] + i[1:] + WT[3:5]
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0:
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0]) 
				epistasis = f_AB - (f_A*f_B)
				e2_3.append(epistasis)

def Ep_e2_4():
	for i in randomizer_list:
		AB = WT[0] + i[:1] + WT[2] + i[1:] + WT[4]     
		if AB == WT:
			e2_4.append(s)
		elif AB[1] == WT[1]:
			e2_4.append(s)
		elif AB[3] == WT[3]:
			e2_4.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe
			f_AB = float(dual_fitness[1].values[0])       # value of fitness corresponding to the sequence found and outputted as float
			A = WT[0] + i[:1] + WT[2:5]
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B = WT[0:3] + i[1:] + WT[4]
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0:
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0])
				epistasis = f_AB - (f_A*f_B)
				e2_4.append(epistasis)

def Ep_e2_5():
	for i in randomizer_list:
		AB = WT[0] + i[:1] + WT[2:4] + i[1:]         
		if AB == WT:
			e2_5.append(s)
		elif AB[1] == WT[1]:
			e2_5.append(s)
		elif AB[4] == WT[4]:
			e2_5.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe
			f_AB = float(dual_fitness[1].values[0])       # value of fitness corresponding to the sequence found and outputted as float
			A = WT[0] + i[:1] + WT[2:5] 
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B = WT[0:4] + i[1:]   
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0:
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0])
				epistasis = f_AB - (f_A*f_B)
				e2_5.append(epistasis)

def Ep_e3_4():
	for i in randomizer_list:
		AB = WT[0:2] + i + WT[4]      
		if AB == WT:
			e3_4.append(s)
		elif AB[2] == WT[2]:
			e3_4.append(s)
		elif AB[3] == WT[3]:
			e3_4.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe
			f_AB = float(dual_fitness[1].values[0])       # value of fitness corresponding to the sequence found and outputted as float
			A = WT[0:2] + i[:1] + WT[3:5] 
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B =  WT[0:3] + i[1:] + WT[4]
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0:
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0])
				epistasis = f_AB - (f_A*f_B)
				e3_4.append(epistasis)

def Ep_e3_5():
	for i in randomizer_list:
		AB = WT[0:2] + i[:1] + WT[3] + i[1:]       
		if AB == WT:
			e3_5.append(s)
		elif AB[2] == WT[2]:
			e3_5.append(s)
		elif AB[4] == WT[4]:
			e3_5.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe
			f_AB = float(dual_fitness[1].values[0])       # value of fitness corresponding to the sequence found and outputted as float
			A = WT[0:2] + i[:1] + WT[3:5]
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B = WT[0:4] + i[1:]  
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0:
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0])
				epistasis = f_AB - (f_A*f_B)
				e3_5.append(epistasis)


def Ep_e4_5():
	for i in randomizer_list:
		AB = WT[0:3] + i      
		if AB == WT:
			e4_5.append(s)
		elif AB[3] == WT[3]:
			e4_5.append(s)
		elif AB[4] == WT[4]:
			e4_5.append(s)
		else:
			dual_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(AB)]      # search for sequence in dataframe
			f_AB = float(dual_fitness[1].values[0])       # value of fitness corresponding to the sequence found and outputted as float
			A =  WT[0:3] + i[:1] + WT[4]
			A_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)] 
			B = WT[0:4] + i[1:]   
			B_fitness = Mutants_fitness_matrix[Mutants_fitness_matrix[0].str.contains(A)]
			if len(A_fitness) > 0 and len(B_fitness) > 0 and len(dual_fitness) > 0:
				f_A = float(A_fitness[1].values[0])        
				f_B = float(B_fitness[1].values[0])
				epistasis = f_AB - (f_A*f_B)
				e4_5.append(epistasis)


Ep_e1_2()
Ep_e1_3()
Ep_e1_4()
Ep_e1_5()
Ep_e2_3()
Ep_e2_4()
Ep_e2_5()
Ep_e3_4()
Ep_e3_5()
Ep_e4_5()

 

# write to a csv file with each of the above clusters being a column

row1 = [s, s, s, s]
row2 = [s, s, s, s]
row3 = [s, s, s, s]
row4 = [s, s, s, s]
row5 = [e1_2[0], e1_2[4], e1_2[8], e1_2[12]]
row6 = [e1_2[1], e1_2[5], e1_2[9], e1_2[13]]
row7 = [e1_2[2], e1_2[6], e1_2[10], e1_2[14]]
row8 = [e1_2[3], e1_2[7], e1_2[11], e1_2[15]]
row9 = [e1_3[0], e1_3[4], e1_3[8], e1_3[12], e2_3[0], e2_3[4], e2_3[8], e2_3[12]]
row10 = [e1_3[1], e1_3[5], e1_3[9], e1_3[13], e2_3[1], e2_3[5], e2_3[9], e2_3[13]]
row11 = [e1_3[2], e1_3[6], e1_3[10], e1_3[14], e2_3[2], e2_3[6], e2_3[10], e2_3[14]]
row12 = [e1_3[3], e1_3[7], e1_3[11], e1_3[15], e2_3[3], e2_3[7], e2_3[11], e2_3[15]]
row13 = [e1_4[0], e1_4[4], e1_4[8], e1_4[12], e2_4[0], e2_4[4], e2_4[8], e2_4[12], e3_4[0], e3_4[4], e3_4[8], e3_4[12]]
row14 = [e1_4[1], e1_4[5], e1_4[9], e1_4[13], e2_4[1], e2_4[5], e2_4[9], e2_4[13], e3_4[1], e3_4[5], e3_4[9], e3_4[13]]
row15 = [e1_4[2], e1_4[6], e1_4[10], e1_4[14], e2_4[2], e2_4[6], e2_4[10], e2_4[14], e3_4[2], e3_4[6], e3_4[10], e3_4[14]]
row16 = [e1_4[3], e1_4[7], e1_4[11], e1_4[15], e2_4[3], e2_4[7], e2_4[11], e2_4[15], e3_4[3], e3_4[7], e3_4[11], e3_4[15]]
row17 = [e1_5[0], e1_5[4], e1_5[8], e1_5[12], e2_5[0], e2_5[4], e2_5[8], e2_5[12], e3_4[0], e3_4[4], e3_4[8], e3_4[12], e4_5[0], e4_5[4], e4_5[8], e4_5[12]]
row18 = [e1_5[1], e1_5[5], e1_5[9], e1_5[13], e2_5[1], e2_5[5], e2_5[9], e2_5[13], e3_4[1], e3_4[5], e3_4[9], e3_4[13], e4_5[1], e4_5[5], e4_5[9], e4_5[13]]
row19 = [e1_5[2], e1_5[6], e1_5[10], e1_5[14], e2_5[2], e2_5[6], e2_5[10], e2_5[14], e3_4[2], e3_4[6], e3_4[10], e3_4[14], e4_5[2], e4_5[6], e4_5[10], e4_5[14]]
row20 = [e1_5[3], e1_5[7], e1_5[11], e1_5[15], e2_5[3], e2_5[7], e2_5[11], e2_5[15], e3_4[3], e3_4[7], e3_4[11], e3_4[15], e4_5[3], e4_5[7], e4_5[11], e4_5[15]]

 

with open('A_epistasis_all.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows([row1])
    writer.writerows([row2])
    writer.writerows([row3])
    writer.writerows([row4])
    writer.writerows([row5])
    writer.writerows([row6])
    writer.writerows([row7])
    writer.writerows([row8])
    writer.writerows([row9])
    writer.writerows([row10])
    writer.writerows([row11])
    writer.writerows([row12])
    writer.writerows([row13])
    writer.writerows([row14])
    writer.writerows([row15])
    writer.writerows([row16])
    writer.writerows([row17])
    writer.writerows([row18])
    writer.writerows([row19])
    writer.writerows([row20])


csvFile.close()



	 






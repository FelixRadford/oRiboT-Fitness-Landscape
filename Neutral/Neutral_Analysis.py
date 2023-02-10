from Bio import SeqIO
from collections import Counter
import csv
import itertools
import pandas as pd
import os

#Obtain the directory where file to be analyzes is located
file_data = open("address.txt", "r")
directory = file_data.read().replace('\n', '')
file_data.close()
 
#Create the name for the file in the correct directory
fName = [f for f in os.listdir(directory) if f.endswith('.fastq')]   # the filename is found within directory and will be used for the next analysis steps
 
filename = os.path.join(directory, fName[0])

good_reads = (rec for rec in \
              SeqIO.parse(filename, "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 20)


primer_reads = (rec for rec in \
                good_reads \
                if "ATGTATCTAAGTAC" in rec)
count = SeqIO.write(primer_reads, "forwardprimerlist.fasta", "fasta")


with open('forwardprimerlist.fasta') as fasta_file:
    sequences = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):   
        Sqnce = str(seq_record.seq)
        sequences.append(Sqnce)

MutRegion = []                                     #all mutant sequences  
for i in sequences:
    start = i.find('TTCATATCGACGGCGGTGTT') + 30   
    end = i.find('TTCATATCGACGGCGGTGTT') + 35 
    part = i[start:end]
    if len(part) == 5:          #Specify the length of the mutagenized region
        MutRegion.append(part)

total_reads = len(MutRegion)

 
WTseq = []
for i in MutRegion:
	if i == 'ATGTC':    #WT sequence of this region
	    WTseq.append(i)

WT_reads = len(WTseq)

MutRegion2 = MutRegion
for i in WTseq:
    MutRegion2.remove(i)
      
unique_reads = len(set(MutRegion))

#Save analysis for this library
text_file = open("Analysis.txt", "w")
text_file.write("total reads: %s\n" % total_reads)
text_file.write("WT reads: %s\n" % WT_reads)
text_file.write("unique reads: %s\n" % unique_reads)
text_file.close()


J = set(MutRegion)
Mt = list(J)


List2 = []
Cntx = Counter(MutRegion)

with open('dict.csv', 'wb') as csv_file:     #find number of repeats for each unique sequence
    writer = csv.writer(csv_file)
    for key, value in Cntx.items():
       writer.writerow([key, value])
       List2.append(value)
Cntx2 = Counter(List2)


with open('dict2.csv', 'wb') as csv_file:    
    writer = csv.writer(csv_file)
    for key, value in Cntx2.items():
       writer.writerow([key, value])
 

with open('MutRegion.csv', 'wb') as csv_file:     
    writer = csv.writer(csv_file)
    for i in MutRegion:
        writer.writerow([i])   


with open('Mutants.csv', 'wb') as csv_file:    #only unique reads 
    writer = csv.writer(csv_file)
    for i in MutRegion2:
        writer.writerow([i])



#generate list of unique 5N sequences:

nucleotide_list = ['A', 'T', 'G', 'C']
N5_sequences = itertools.product(nucleotide_list, repeat = 5)
N5_sequences = [''.join(i) for i in itertools.product(nucleotide_list, repeat = 5)]

N5_all = pd.DataFrame.from_dict(N5_sequences)
N5_all.columns = ['index']
Read_freqs_data =  pd.DataFrame.from_dict(Cntx, orient='index')
Read_freqs_data.reset_index(level=0, inplace=True)


N5_calc = pd.merge(N5_all, Read_freqs_data, on='index', how='left')    
N5_intersection = pd.merge(N5_all, Read_freqs_data, on='index')    

N5_intersection.to_csv('N5_dict.csv')

#calculate percentages of total for each unique sequence
summ = N5_calc[0].sum()
percentages = (N5_calc[0]/summ)*100
N5_calc['percentages']= percentages


#write dataframe into file, containing all unique sequences in order and counts for each one
N5_calc.to_csv('unique_sequences.csv')
N5_intersection.to_csv('intersect_sequences.csv')

A1 = []
C1 = []
T1 = []
G1 = []
A2 = []
C2 = []
T2 = []
G2 = []
A3 = []
C3 = []
T3 = []
G3 = []
A4 = []
C4 = []
T4 = []
G4 = []
A5 = []
C5 = []
T5 = []
G5 = []
A6 = []
C6 = []
T6 = []
G6 = []
A7 = []
C7 = []
T7 = []
G7 = []

def Posn_1():
    for i in MutRegion2[1:]:
        if i[0] == 'A':
            A1.append(1)
        if i[0] == 'C':
            C1.append(1)
        if i[0] == 'T':
            T1.append(1)
        if i[0] == 'G':
            G1.append(1)
        else:
            pass

def Posn_2():
    for i in MutRegion2[1:]:
        if i[1] == 'A':
            A2.append(1)
        if i[1] == 'C':
            C2.append(1)
        if i[1] == 'T':
            T2.append(1)
        if i[1] == 'G':
            G2.append(1)
        else:
            pass
    
def Posn_3():
    for i in MutRegion2[1:]:
        if i[2] == 'A':
            A3.append(1)
        if i[2] == 'C':
            C3.append(1)
        if i[2] == 'T':
            T3.append(1)
        if i[2] == 'G':
            G3.append(1)
        else:
            pass
    
def Posn_4():
    for i in MutRegion2[1:]:
        if i[3] == 'A':
            A4.append(1)
        if i[3] == 'C':
            C4.append(1)
        if i[3] == 'T':
            T4.append(1)
        if i[3] == 'G':
            G4.append(1)
        else:
            pass


def Posn_5():
    for i in MutRegion2[1:]:
        if i[4] == 'A':
            A5.append(1)
        if i[4] == 'C':
            C5.append(1)
        if i[4] == 'T':
            T5.append(1)
        if i[4] == 'G':
            G5.append(1)
        else:
            pass 

def Posn_6():
    for i in MutRegion2[1:]:
        if i[5] == 'A':
            A6.append(1)
        if i[5] == 'C':
            C6.append(1)
        if i[5] == 'T':
            T6.append(1)
        if i[5] == 'G':
            G6.append(1)
        else:
            pass 

def Posn_7():
    for i in MutRegion2[1:]:
        if i[6] == 'A':
            A7.append(1)
        if i[6] == 'C':
            C7.append(1)
        if i[6] == 'T':
            T7.append(1)
        if i[6] == 'G':
            G7.append(1)
        else:
            pass 
Posn_1()
Posn_2()
Posn_3()
Posn_4()
Posn_5()


A1 = len(A1)
C1 = len(C1)
T1 = len(T1)
G1 = len(G1)
A2 = len(A2)
C2 = len(C2)
T2 = len(T2)
G2 = len(G2)
A3 = len(A3)
C3 = len(C3)
T3 = len(T3)
G3 = len(G3)
A4 = len(A4)
C4 = len(C4)
T4 = len(T4)
G4 = len(G4)
A5 = len(A5)
C5 = len(C5)
T5 = len(T5)
G5 = len(G5)
A6 = len(A6)
C6 = len(C6)
T6 = len(T6)
G6 = len(G6)
A7 = len(A7)
C7 = len(C7)
T7 = len(T7)
G7 = len(G7)

row1 = [A1, A2, A3, A4, A5]
row2 = [C1, C2, C3, C4, C5]
row3 = [T1, T2, T3, T4, T5]
row4 = [G1, G2, G3, G4, G5]


with open('Sites-Matrix.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows([row1])
    writer.writerows([row2])
    writer.writerows([row3])
    writer.writerows([row4])


csvFile.close()

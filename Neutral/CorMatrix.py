import csv
import matplotlib.pyplot as plt
import pylab
import pandas

#Set up lists for correlation matrix counts------------------------------------

A_A_1_2 = []
A_C_1_2 = []
A_T_1_2 = []
A_G_1_2 = []
C_C_1_2 = []
C_T_1_2 = []
C_G_1_2 = []
T_T_1_2 = []
T_G_1_2 = []
G_G_1_2 = []
C_A_1_2 = []
T_A_1_2 = []
G_A_1_2 = []
T_C_1_2 = []
G_C_1_2 = []
G_T_1_2 = []

A_A_1_3 = []
A_C_1_3 = []
A_T_1_3 = []
A_G_1_3 = []
C_C_1_3 = []
C_T_1_3 = []
C_G_1_3 = []
T_T_1_3 = []
T_G_1_3 = []
G_G_1_3 = []
C_A_1_3 = []
T_A_1_3 = []
G_A_1_3 = []
T_C_1_3 = []
G_C_1_3 = []
G_T_1_3 = []

A_A_1_4 = []
A_C_1_4 = []
A_T_1_4 = []
A_G_1_4 = []
C_C_1_4 = []
C_T_1_4 = []
C_G_1_4 = []
T_T_1_4 = []
T_G_1_4 = []
G_G_1_4 = []
C_A_1_4 = []
T_A_1_4 = []
G_A_1_4 = []
T_C_1_4 = []
G_C_1_4 = []
G_T_1_4 = []

A_A_1_5 = []
A_C_1_5 = []
A_T_1_5 = []
A_G_1_5 = []
C_C_1_5 = []
C_T_1_5 = []
C_G_1_5 = []
T_T_1_5 = []
T_G_1_5 = []
G_G_1_5 = []
C_A_1_5 = []
T_A_1_5 = []
G_A_1_5 = []
T_C_1_5 = []
G_C_1_5 = []
G_T_1_5 = []

A_A_2_3 = []
A_C_2_3 = []
A_T_2_3 = []
A_G_2_3 = []
C_C_2_3 = []
C_T_2_3 = []
C_G_2_3 = []
T_T_2_3 = []
T_G_2_3 = []
G_G_2_3 = []
C_A_2_3 = []
T_A_2_3 = []
G_A_2_3 = []
T_C_2_3 = []
G_C_2_3 = []
G_T_2_3 = []

A_A_2_4 = []
A_C_2_4 = []
A_T_2_4 = []
A_G_2_4 = []
C_C_2_4 = []
C_T_2_4 = []
C_G_2_4 = []
T_T_2_4 = []
T_G_2_4 = []
G_G_2_4 = []
C_A_2_4 = []
T_A_2_4 = []
G_A_2_4 = []
T_C_2_4 = []
G_C_2_4 = []
G_T_2_4 = []

A_A_2_5 = []
A_C_2_5 = []
A_T_2_5 = []
A_G_2_5 = []
C_C_2_5 = []
C_T_2_5 = []
C_G_2_5 = []
T_T_2_5 = []
T_G_2_5 = []
G_G_2_5 = []
C_A_2_5 = []
T_A_2_5 = []
G_A_2_5 = []
T_C_2_5 = []
G_C_2_5 = []
G_T_2_5 = []

A_A_3_4 = []
A_C_3_4 = []
A_T_3_4 = []
A_G_3_4 = []
C_C_3_4 = []
C_T_3_4 = []
C_G_3_4 = []
T_T_3_4 = []
T_G_3_4 = []
G_G_3_4 = []
C_A_3_4 = []
T_A_3_4 = []
G_A_3_4 = []
T_C_3_4 = []
G_C_3_4 = []
G_T_3_4 = []

A_A_3_5 = []
A_C_3_5 = []
A_T_3_5 = []
A_G_3_5 = []
C_C_3_5 = []
C_T_3_5 = []
C_G_3_5 = []
T_T_3_5 = []
T_G_3_5 = []
G_G_3_5 = []
C_A_3_5 = []
T_A_3_5 = []
G_A_3_5 = []
T_C_3_5 = []
G_C_3_5 = []
G_T_3_5 = []

A_A_4_5 = []
A_C_4_5 = []
A_T_4_5 = []
A_G_4_5 = []
C_C_4_5 = []
C_T_4_5 = []
C_G_4_5 = []
T_T_4_5 = []
T_G_4_5 = []
G_G_4_5 = []
C_A_4_5 = []
T_A_4_5 = []
G_A_4_5 = []
T_C_4_5 = []
G_C_4_5 = []
G_T_4_5 = []

#Functions to split the strings and populate the above lists------------------------------------

#read the total mutants (MutRegion2 written to a file, Mutants.csv); populate into list called Mutants

df = pandas.read_csv('Mutants.csv', header=None)
Mutants = df[0].tolist()


def Split1_2():  #split strings of 5N into 1st and 2nd joined nts, and populated into appropriate lists
    for i in Mutants:
        part = i[0:2]
        if part == 'AA':
            A_A_1_2.append(part)
        if part == 'AC':
        	A_C_1_2.append(part)
        if part == 'AT':
        	A_T_1_2.append(part)
        if part == 'AG':
        	A_G_1_2.append(part)
        if part == 'CC':
        	C_C_1_2.append(part)
        if part == 'CT':
        	C_T_1_2.append(part)
        if part == 'CG':
        	C_G_1_2.append(part)
        if part == 'TT':
        	T_T_1_2.append(part)
        if part == 'TG':
        	T_G_1_2.append(part)
        if part == 'GG':
        	G_G_1_2.append(part)
        if part == 'CA':
        	C_A_1_2.append(part)
        if part == 'TA':
        	T_A_1_2.append(part)
        if part == 'GA':
        	G_A_1_2.append(part)
        if part == 'TC':
        	T_C_1_2.append(part)
        if part == 'GC':
        	G_C_1_2.append(part)
        if part == 'GT':
        	G_T_1_2.append(part)

def Split1_3():  
    for i in Mutants:
        part = i[0] + i[2]
        if part == 'AA':
            A_A_1_3.append(part)
        if part == 'AC':
        	A_C_1_3.append(part)
        if part == 'AT':
        	A_T_1_3.append(part)
        if part == 'AG':
        	A_G_1_3.append(part)
        if part == 'CC':
        	C_C_1_3.append(part)
        if part == 'CT':
        	C_T_1_3.append(part)
        if part == 'CG':
        	C_G_1_3.append(part)
        if part == 'TT':
        	T_T_1_3.append(part)
        if part == 'TG':
        	T_G_1_3.append(part)
        if part == 'GG':
        	G_G_1_3.append(part)
        if part == 'CA':
        	C_A_1_3.append(part)
        if part == 'TA':
        	T_A_1_3.append(part)
        if part == 'GA':
        	G_A_1_3.append(part)
        if part == 'TC':
        	T_C_1_3.append(part)
        if part == 'GC':
        	G_C_1_3.append(part)
        if part == 'GT':
        	G_T_1_3.append(part)


def Split1_4():  
    for i in Mutants:
        part = i[0] + i[3]
        if part == 'AA':
            A_A_1_4.append(part)
        if part == 'AC':
        	A_C_1_4.append(part)
        if part == 'AT':
        	A_T_1_4.append(part)
        if part == 'AG':
        	A_G_1_4.append(part)
        if part == 'CC':
        	C_C_1_4.append(part)
        if part == 'CT':
        	C_T_1_4.append(part)
        if part == 'CG':
        	C_G_1_4.append(part)
        if part == 'TT':
        	T_T_1_4.append(part)
        if part == 'TG':
        	T_G_1_4.append(part)
        if part == 'GG':
        	G_G_1_4.append(part)
        if part == 'CA':
        	C_A_1_4.append(part)
        if part == 'TA':
        	T_A_1_4.append(part)
        if part == 'GA':
        	G_A_1_4.append(part)
        if part == 'TC':
        	T_C_1_4.append(part)
        if part == 'GC':
        	G_C_1_4.append(part)
        if part == 'GT':
        	G_T_1_4.append(part)

def Split1_5():  
    for i in Mutants:
        part = i[0] + i[4]
        if part == 'AA':
            A_A_1_5.append(part)
        if part == 'AC':
        	A_C_1_5.append(part)
        if part == 'AT':
        	A_T_1_5.append(part)
        if part == 'AG':
        	A_G_1_5.append(part)
        if part == 'CC':
        	C_C_1_5.append(part)
        if part == 'CT':
        	C_T_1_5.append(part)
        if part == 'CG':
        	C_G_1_5.append(part)
        if part == 'TT':
        	T_T_1_5.append(part)
        if part == 'TG':
        	T_G_1_5.append(part)
        if part == 'GG':
        	G_G_1_5.append(part)
        if part == 'CA':
        	C_A_1_5.append(part)
        if part == 'TA':
        	T_A_1_5.append(part)
        if part == 'GA':
        	G_A_1_5.append(part)
        if part == 'TC':
        	T_C_1_5.append(part)
        if part == 'GC':
        	G_C_1_5.append(part)
        if part == 'GT':
        	G_T_1_5.append(part)

def Split2_3():  
    for i in Mutants:
        part = i[1] + i[2]
        if part == 'AA':
            A_A_2_3.append(part)
        if part == 'AC':
        	A_C_2_3.append(part)
        if part == 'AT':
        	A_T_2_3.append(part)
        if part == 'AG':
        	A_G_2_3.append(part)
        if part == 'CC':
        	C_C_2_3.append(part)
        if part == 'CT':
        	C_T_2_3.append(part)
        if part == 'CG':
        	C_G_2_3.append(part)
        if part == 'TT':
        	T_T_2_3.append(part)
        if part == 'TG':
        	T_G_2_3.append(part)
        if part == 'GG':
        	G_G_2_3.append(part)
        if part == 'CA':
        	C_A_2_3.append(part)
        if part == 'TA':
        	T_A_2_3.append(part)
        if part == 'GA':
        	G_A_2_3.append(part)
        if part == 'TC':
        	T_C_2_3.append(part)
        if part == 'GC':
        	G_C_2_3.append(part)
        if part == 'GT':
        	G_T_2_3.append(part)

def Split2_4():  
    for i in Mutants:
        part = i[1] + i[3]
        if part == 'AA':
            A_A_2_4.append(part)
        if part == 'AC':
        	A_C_2_4.append(part)
        if part == 'AT':
        	A_T_2_4.append(part)
        if part == 'AG':
        	A_G_2_4.append(part)
        if part == 'CC':
        	C_C_2_4.append(part)
        if part == 'CT':
        	C_T_2_4.append(part)
        if part == 'CG':
        	C_G_2_4.append(part)
        if part == 'TT':
        	T_T_2_4.append(part)
        if part == 'TG':
        	T_G_2_4.append(part)
        if part == 'GG':
        	G_G_2_4.append(part)
        if part == 'CA':
        	C_A_2_4.append(part)
        if part == 'TA':
        	T_A_2_4.append(part)
        if part == 'GA':
        	G_A_2_4.append(part)
        if part == 'TC':
        	T_C_2_4.append(part)
        if part == 'GC':
        	G_C_2_4.append(part)
        if part == 'GT':
        	G_T_2_4.append(part)


def Split2_5():  
    for i in Mutants:
        part = i[1] + i[4]
        if part == 'AA':
            A_A_2_5.append(part)
        if part == 'AC':
        	A_C_2_5.append(part)
        if part == 'AT':
        	A_T_2_5.append(part)
        if part == 'AG':
        	A_G_2_5.append(part)
        if part == 'CC':
        	C_C_2_5.append(part)
        if part == 'CT':
        	C_T_2_5.append(part)
        if part == 'CG':
        	C_G_2_5.append(part)
        if part == 'TT':
        	T_T_2_5.append(part)
        if part == 'TG':
        	T_G_2_5.append(part)
        if part == 'GG':
        	G_G_2_5.append(part)
        if part == 'CA':
        	C_A_2_5.append(part)
        if part == 'TA':
        	T_A_2_5.append(part)
        if part == 'GA':
        	G_A_2_5.append(part)
        if part == 'TC':
        	T_C_2_5.append(part)
        if part == 'GC':
        	G_C_2_5.append(part)
        if part == 'GT':
        	G_T_2_5.append(part)

def Split3_4():  
    for i in Mutants:
        part = i[2] + i[3]
        if part == 'AA':
            A_A_3_4.append(part)
        if part == 'AC':
        	A_C_3_4.append(part)
        if part == 'AT':
        	A_T_3_4.append(part)
        if part == 'AG':
        	A_G_3_4.append(part)
        if part == 'CC':
        	C_C_3_4.append(part)
        if part == 'CT':
        	C_T_3_4.append(part)
        if part == 'CG':
        	C_G_3_4.append(part)
        if part == 'TT':
        	T_T_3_4.append(part)
        if part == 'TG':
        	T_G_3_4.append(part)
        if part == 'GG':
        	G_G_3_4.append(part)
        if part == 'CA':
        	C_A_3_4.append(part)
        if part == 'TA':
        	T_A_3_4.append(part)
        if part == 'GA':
        	G_A_3_4.append(part)
        if part == 'TC':
        	T_C_3_4.append(part)
        if part == 'GC':
        	G_C_3_4.append(part)
        if part == 'GT':
        	G_T_3_4.append(part)


def Split3_5():  
    for i in Mutants:
        part = i[2] + i[4]
        if part == 'AA':
            A_A_3_5.append(part)
        if part == 'AC':
        	A_C_3_5.append(part)
        if part == 'AT':
        	A_T_3_5.append(part)
        if part == 'AG':
        	A_G_3_5.append(part)
        if part == 'CC':
        	C_C_3_5.append(part)
        if part == 'CT':
        	C_T_3_5.append(part)
        if part == 'CG':
        	C_G_3_5.append(part)
        if part == 'TT':
        	T_T_3_5.append(part)
        if part == 'TG':
        	T_G_3_5.append(part)
        if part == 'GG':
        	G_G_3_5.append(part)
        if part == 'CA':
        	C_A_3_5.append(part)
        if part == 'TA':
        	T_A_3_5.append(part)
        if part == 'GA':
        	G_A_3_5.append(part)
        if part == 'TC':
        	T_C_3_5.append(part)
        if part == 'GC':
        	G_C_3_5.append(part)
        if part == 'GT':
        	G_T_3_5.append(part)


def Split4_5():  
    for i in Mutants:
        part = i[3] + i[4]
        if part == 'AA':
            A_A_4_5.append(part)
        if part == 'AC':
        	A_C_4_5.append(part)
        if part == 'AT':
        	A_T_4_5.append(part)
        if part == 'AG':
        	A_G_4_5.append(part)
        if part == 'CC':
        	C_C_4_5.append(part)
        if part == 'CT':
        	C_T_4_5.append(part)
        if part == 'CG':
        	C_G_4_5.append(part)
        if part == 'TT':
        	T_T_4_5.append(part)
        if part == 'TG':
        	T_G_4_5.append(part)
        if part == 'GG':
        	G_G_4_5.append(part)
        if part == 'CA':
        	C_A_4_5.append(part)
        if part == 'TA':
        	T_A_4_5.append(part)
        if part == 'GA':
        	G_A_4_5.append(part)
        if part == 'TC':
        	T_C_4_5.append(part)
        if part == 'GC':
        	G_C_4_5.append(part)
        if part == 'GT':
        	G_T_4_5.append(part)

# Execute the functions

Split1_2()
Split1_3()
Split1_4()
Split1_5()
Split2_3()
Split2_4()
Split2_5()
Split3_4()
Split3_5()
Split4_5()


# compute the lengths of the populated lists

A_A_1_2 = len(A_A_1_2)
A_C_1_2 = len(A_C_1_2)
A_T_1_2 = len(A_T_1_2)
A_G_1_2 = len(A_G_1_2)
C_C_1_2 = len(C_C_1_2)
C_T_1_2 = len(C_T_1_2)
C_G_1_2 = len(C_G_1_2)
T_T_1_2 = len(T_T_1_2)
T_G_1_2 = len(T_G_1_2)
G_G_1_2 = len(G_G_1_2)
C_A_1_2 = len(C_A_1_2)
T_A_1_2 = len(T_A_1_2)
G_A_1_2 = len(G_A_1_2)
T_C_1_2 = len(T_C_1_2)
G_C_1_2 = len(G_C_1_2)
G_T_1_2 = len(G_T_1_2)

A_A_1_3 = len(A_A_1_3)
A_C_1_3 = len(A_C_1_3)
A_T_1_3 = len(A_T_1_3)
A_G_1_3 = len(A_G_1_3)
C_C_1_3 = len(C_C_1_3)
C_T_1_3 = len(C_T_1_3)
C_G_1_3 = len(C_G_1_3)
T_T_1_3 = len(T_T_1_3)
T_G_1_3 = len(T_G_1_3)
G_G_1_3 = len(G_G_1_3)
C_A_1_3 = len(C_A_1_3)
T_A_1_3 = len(T_A_1_3)
G_A_1_3 = len(G_A_1_3)
T_C_1_3 = len(T_C_1_3)
G_C_1_3 = len(G_C_1_3)
G_T_1_3 = len(G_T_1_3)

A_A_1_4 = len(A_A_1_4)
A_C_1_4 = len(A_C_1_4)
A_T_1_4 = len(A_T_1_4)
A_G_1_4 = len(A_G_1_4)
C_C_1_4 = len(C_C_1_4)
C_T_1_4 = len(C_T_1_4)
C_G_1_4 = len(C_G_1_4)
T_T_1_4 = len(T_T_1_4)
T_G_1_4 = len(T_G_1_4)
G_G_1_4 = len(G_G_1_4)
C_A_1_4 = len(C_A_1_4)
T_A_1_4 = len(T_A_1_4)
G_A_1_4 = len(G_A_1_4)
T_C_1_4 = len(T_C_1_4)
G_C_1_4 = len(G_C_1_4)
G_T_1_4 = len(G_T_1_4)

A_A_1_5 = len(A_A_1_5)
A_C_1_5 = len(A_C_1_5)
A_T_1_5 = len(A_T_1_5)
A_G_1_5 = len(A_G_1_5)
C_C_1_5 = len(C_C_1_5)
C_T_1_5 = len(C_T_1_5)
C_G_1_5 = len(C_G_1_5)
T_T_1_5 = len(T_T_1_5)
T_G_1_5 = len(T_G_1_5)
G_G_1_5 = len(G_G_1_5)
C_A_1_5 = len(C_A_1_5)
T_A_1_5 = len(T_A_1_5)
G_A_1_5 = len(G_A_1_5)
T_C_1_5 = len(T_C_1_5)
G_C_1_5 = len(G_C_1_5)
G_T_1_5 = len(G_T_1_5)

A_A_2_3 = len(A_A_2_3)
A_C_2_3 = len(A_C_2_3)
A_T_2_3 = len(A_T_2_3)
A_G_2_3 = len(A_G_2_3)
C_C_2_3 = len(C_C_2_3)
C_T_2_3 = len(C_T_2_3)
C_G_2_3 = len(C_G_2_3)
T_T_2_3 = len(T_T_2_3)
T_G_2_3 = len(T_G_2_3)
G_G_2_3 = len(G_G_2_3)
C_A_2_3 = len(C_A_2_3)
T_A_2_3 = len(T_A_2_3)
G_A_2_3 = len(G_A_2_3)
T_C_2_3 = len(T_C_2_3)
G_C_2_3 = len(G_C_2_3)
G_T_2_3 = len(G_T_2_3)

A_A_2_4 = len(A_A_2_4)
A_C_2_4 = len(A_C_2_4)
A_T_2_4 = len(A_T_2_4)
A_G_2_4 = len(A_G_2_4)
C_C_2_4 = len(C_C_2_4)
C_T_2_4 = len(C_T_2_4)
C_G_2_4 = len(C_G_2_4)
T_T_2_4 = len(T_T_2_4)
T_G_2_4 = len(T_G_2_4)
G_G_2_4 = len(G_G_2_4)
C_A_2_4 = len(C_A_2_4)
T_A_2_4 = len(T_A_2_4)
G_A_2_4 = len(G_A_2_4)
T_C_2_4 = len(T_C_2_4)
G_C_2_4 = len(G_C_2_4)
G_T_2_4 = len(G_T_2_4)

A_A_2_5 = len(A_A_2_5)
A_C_2_5 = len(A_C_2_5)
A_T_2_5 = len(A_T_2_5)
A_G_2_5 = len(A_G_2_5)
C_C_2_5 = len(C_C_2_5)
C_T_2_5 = len(C_T_2_5)
C_G_2_5 = len(C_G_2_5)
T_T_2_5 = len(T_T_2_5)
T_G_2_5 = len(T_G_2_5)
G_G_2_5 = len(G_G_2_5)
C_A_2_5 = len(C_A_2_5)
T_A_2_5 = len(T_A_2_5)
G_A_2_5 = len(G_A_2_5)
T_C_2_5 = len(T_C_2_5)
G_C_2_5 = len(G_C_2_5)
G_T_2_5 = len(G_T_2_5)

A_A_3_4 = len(A_A_3_4)
A_C_3_4 = len(A_C_3_4)
A_T_3_4 = len(A_T_3_4)
A_G_3_4 = len(A_G_3_4)
C_C_3_4 = len(C_C_3_4)
C_T_3_4 = len(C_T_3_4)
C_G_3_4 = len(C_G_3_4)
T_T_3_4 = len(T_T_3_4)
T_G_3_4 = len(T_G_3_4)
G_G_3_4 = len(G_G_3_4)
C_A_3_4 = len(C_A_3_4)
T_A_3_4 = len(T_A_3_4)
G_A_3_4 = len(G_A_3_4)
T_C_3_4 = len(T_C_3_4)
G_C_3_4 = len(G_C_3_4)
G_T_3_4 = len(G_T_3_4)

A_A_3_5 = len(A_A_3_5)
A_C_3_5 = len(A_C_3_5)
A_T_3_5 = len(A_T_3_5)
A_G_3_5 = len(A_G_3_5)
C_C_3_5 = len(C_C_3_5)
C_T_3_5 = len(C_T_3_5)
C_G_3_5 = len(C_G_3_5)
T_T_3_5 = len(T_T_3_5)
T_G_3_5 = len(T_G_3_5)
G_G_3_5 = len(G_G_3_5)
C_A_3_5 = len(C_A_3_5)
T_A_3_5 = len(T_A_3_5)
G_A_3_5 = len(G_A_3_5)
T_C_3_5 = len(T_C_3_5)
G_C_3_5 = len(G_C_3_5)
G_T_3_5 = len(G_T_3_5)

A_A_4_5 = len(A_A_4_5)
A_C_4_5 = len(A_C_4_5)
A_T_4_5 = len(A_T_4_5)
A_G_4_5 = len(A_G_4_5)
C_C_4_5 = len(C_C_4_5)
C_T_4_5 = len(C_T_4_5)
C_G_4_5 = len(C_G_4_5)
T_T_4_5 = len(T_T_4_5)
T_G_4_5 = len(T_G_4_5)
G_G_4_5 = len(G_G_4_5)
C_A_4_5 = len(C_A_4_5)
T_A_4_5 = len(T_A_4_5)
G_A_4_5 = len(G_A_4_5)
T_C_4_5 = len(T_C_4_5)
G_C_4_5 = len(G_C_4_5)
G_T_4_5 = len(G_T_4_5)




# write to a csv file with each of the above clusters being a column

row1 = [A_A_1_2, A_A_1_3, A_A_1_4, A_A_1_5, A_A_2_3, A_A_2_4, A_A_2_5, A_A_3_4, A_A_3_5, A_A_4_5]
row2 = [A_C_1_2, A_C_1_3, A_C_1_4, A_C_1_5, A_C_2_3, A_C_2_4, A_C_2_5, A_C_3_4, A_C_3_5, A_C_4_5]
row3 = [A_T_1_2, A_T_1_3, A_T_1_4, A_T_1_5, A_T_2_3, A_T_2_4, A_T_2_5, A_T_3_4, A_T_3_5, A_T_4_5]
row4 = [A_G_1_2, A_G_1_3, A_G_1_4, A_G_1_5, A_G_2_3, A_G_2_4, A_G_2_5, A_G_3_4, A_G_3_5, A_G_4_5]
row5 = [C_C_1_2, C_C_1_3, C_C_1_4, C_C_1_5, C_C_2_3, C_C_2_4, C_C_2_5, C_C_3_4, C_C_3_5, C_C_4_5]
row6 = [C_T_1_2, C_T_1_3, C_T_1_4, C_T_1_5, C_T_2_3, C_T_2_4, C_T_2_5, C_T_3_4, C_T_3_5, C_T_4_5]
row7 = [C_G_1_2, C_G_1_3, C_G_1_4, C_G_1_5, C_G_2_3, C_G_2_4, C_G_2_5, C_G_3_4, C_G_3_5, C_G_4_5]
row8 = [T_T_1_2, T_T_1_3, T_T_1_4, T_T_1_5, T_T_2_3, T_T_2_4, T_T_2_5, T_T_3_4, T_T_3_5, T_T_4_5]
row9 = [T_G_1_2, T_G_1_3, T_G_1_4, T_G_1_5, T_G_2_3, T_G_2_4, T_G_2_5, T_G_3_4, T_G_3_5, T_G_4_5]
row10 = [G_G_1_2, G_G_1_3, G_G_1_4, G_G_1_5, G_G_2_3, G_G_2_4, G_G_2_5, G_G_3_4, G_G_3_5, G_G_4_5]
row11 = [C_A_1_2, C_A_1_3, C_A_1_4, C_A_1_5, C_A_2_3, C_A_2_4, C_A_2_5, C_A_3_4, A_C_3_5, C_A_4_5]      
row12 = [T_A_1_2, T_A_1_3, T_A_1_4, T_A_1_5, T_A_2_3, T_A_2_4, T_A_2_5, T_A_3_4, T_A_3_5, T_A_4_5]
row13 = [G_A_1_2, G_A_1_3, G_A_1_4, G_A_1_5, G_A_2_3, G_A_2_4, G_A_2_5, G_A_3_4, G_A_3_5, G_A_4_5]
row14 = [T_C_1_2, T_C_1_3, T_C_1_4, T_C_1_5, T_C_2_3, T_C_2_4, T_C_2_5, T_C_3_4, T_C_3_5, T_C_4_5]
row15 = [G_C_1_2, G_C_1_3, G_C_1_4, G_C_1_5, G_C_2_3, G_C_2_4, G_C_2_5, G_C_3_4, G_C_3_5, G_C_4_5]
row16 = [G_T_1_2, G_T_1_3, G_T_1_4, G_T_1_5, G_T_2_3, G_T_2_4, G_T_2_5, G_T_3_4, G_T_3_5, G_T_4_5]


with open('cormatrix.csv', 'wb') as csvFile:
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

csvFile.close()

























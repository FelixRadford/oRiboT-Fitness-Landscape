# oRiboT-Fitness-Landscape

Python scripts used for the publication 
**Mapping the in vivo fitness landscape of a tethered ribosome** 
by Felix Radford, Jess Rinehart, and Farren J. Isaacs.


## Table of Contents
- [Prerequisites](#prerequisites)
- [Usage](#usage)




## Prerequisites 

- [Python 3.7](https://www.python.org/)  
- [scipy](https://anaconda.org/anaconda/scipy) package for Python
- [pandas](https://anaconda.org/anaconda/pandas) package for Python 
- [numpy](https://anaconda.org/anaconda/pandas) package for Python 
- [seaborn](https://anaconda.org/anaconda/seaborn) package for Python 
- [matplotlib](https://anaconda.org/anaconda/matplotlib) package for Python 
- [biopython](https://anaconda.org/anaconda/biopython) package for Python 


## Usage

A .fastq file containing the pooled positive and negative selections, and unselected population, can be analyzed by running Gen_Analysis.py. This will generate the analysis of each library, as well as an enrichment matrix for all sites of the library, in the folders labeled “Negative”, “Neutral”, and “Positive.” The example script here is shown for the 2503-2507 library. The WT sequence of the library should be set to the corresponding value for each PTC library analyzed within Negative_Analysis.py, Positive_Analysis.py, and Neutral_Analysis.py. For example, it is 'ATGTC' for the 2503-2507 library. 

To generate covariation matrices of all positions in a library, run Gen_Analysis_Cormatrix.py, which will create .csv files of the covariation matrices in “Negative”, “Neutral”, and “Positive” folders, as well as a finalized covariation matrix for the library.

Fitness of each ribosome mutant was calculated by the following equation: f(A) = (Ap – Ao)/(An – Ao), where where f(A) is the fitness of each oRiboT mutant, and the enrichment of each oRiboT mutant (determined by quantification of the NGS reads as detailed above) is given by Ap in the positive selection, An in the negative selection, and A0 in the not selected population. Running fitness_calc.py will calculate the fitness of each oRiboT mutant in the population. 

We quantified the flexibility to mutation at each position in the population by computing the [Shannon entropy](https://onlinelibrary.wiley.com/doi/10.1002/j.1538-7305.1948.tb01338.x). Running entropy_all_positions.py calculates the entropy at each position in the library. 

We calculated epistasis for each PTC library using epistasis_matrix.py. Epistasis was calculated based on a [previous work](https://pubmed.ncbi.nlm.nih.gov/27080104/)analyzing the epistasis between the bases of an RNA element according to the following equation: ε = f(ABi) – f(Ai)*f(Bi), where ε is the epistasis, f(A) is the fitness of point mutant with point mutation A, f(B) is the fitness of point mutant with point mutation B, and f(AB) is the fitness of a mutant having both mutations A and B.

## Authors

* **Felix Radford** 

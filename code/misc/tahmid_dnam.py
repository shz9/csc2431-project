import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# load data
# data = pd.read_table("~/data/charge_kabuki_matrix.txt", header=0, index_col=0)
# load annotations for each patient
patient_annotations = pd.read_table("C:/Users/Moses Lab/Documents/Data/DNAm_csc2431/patient_annotations.txt",
                                    index_col=0)
diseases = patient_annotations.iloc[3]
types = patient_annotations.iloc[2]
# genomic coordinates of CpGs

print patient_annotations.describe()
print "----"
print patient_annotations.columns
print "----"
print diseases
print types
"""coords = pd.read_table("data/id2coord.txt", header=0, index_col=0)
# remove X & Y chromosomes
sexChr = coords[(coords['Chromosome_36'] == "Y") | (coords['Chromosome_36'] == "X")].index
data = data.iloc[[(i not in sexChr) for i in data.index]]

X = data.transpose()
# remove CpGs with missing betas
X = X.dropna(axis=1, how='any')"""

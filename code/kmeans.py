import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pickle
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import LabelEncoder
from sklearn.semi_supervised import LabelSpreading
import random


def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


def filter_probes(data, probes):
    return data[[probe for probe in data.columns if probe in probes]]


def kmeansWrapper(data, y, k, cores=1):
    kmeans = KMeans(n_clusters=k, n_jobs=cores).fit(data)
    ari = adjusted_rand_score(y, kmeans.labels_)
    return ari


def LabelSpreadingWrapper(Xtrain, ytrain, Xtest, unlabeled):
    s = random.sample(range(len(ytrain)), int(unlabeled))
    ytrain[s] = -1

    lp = LabelSpreading(kernel='knn', n_neighbors=7, alpha=0.2)
    lp.fit(Xtrain, ytrain)
    probs = lp.predict_proba(Xtest)
    return probs


with open('ck.pkl', 'rb') as f:
    X = pickle.load(f)

Xdisc = X.iloc[[i for i in range(len(X)) if 'CHD7 LOF' in X.index[i] or 'KMT2D LOF' in X.index[i] or 110 <= i < 150]]
diseases = [find_between(sample, 'blood ', ' genomic') for sample in Xdisc.index]
le = LabelEncoder()
labels = le.fit_transform(diseases)
# differentially methylated
chd7Sigs = pd.read_csv("R/chd7sigs.csv")
kmt2dSigs = pd.read_csv("R/kmt2dsigs.csv")
chd7Sigs = chd7Sigs['TargetID'].values
kmt2dSigs = kmt2dSigs['TargetID'].values
sigs = np.unique(np.append(chd7Sigs, kmt2dSigs))
XsigsDisc = filter_probes(Xdisc, sigs)

# statistical testing
ttest = pd.read_csv("t_test_results.csv")
utest = pd.read_csv("u_test_results.csv")
ttestProbes = ttest['Probe_ID'].values
utestProbes = utest['Probe_ID'].values
unionProbes = np.unique(np.append(ttestProbes, utestProbes))
interProbes = np.intersect1d(ttestProbes, utestProbes)
XttestDisc = filter_probes(Xdisc, ttestProbes)
XutestDisc = filter_probes(Xdisc, utestProbes)
XunionDisc = filter_probes(Xdisc, unionProbes)
XinterDisc = filter_probes(Xdisc, interProbes)

# PCA | mean: 0.27922529309066607, std: 0.016856514055956354
pca = PCA(n_components=40, random_state=0)
Xpca = pca.fit_transform(X)
Xpca = pd.DataFrame(Xpca)
Xpca.index = X.index
XpcaDisc = Xpca.iloc[[i for i in range(len(Xpca)) if Xpca.index[i] in Xdisc.index]]

n = 100
ari = [0]*n
for i in range(n):
    ari[i] = kmeansWrapper(XpcaDisc, labels, 3, cores=3)


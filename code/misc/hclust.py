import pandas as pd
import numpy as np
import pickle
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.preprocessing import LabelEncoder
import seaborn as sns


def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


def filter_probes(data, probes):
    return data[[probe for probe in data.columns if probe in probes]]


def hclustWrapper(data, y, k, affinity, linkage):
    ac = AgglomerativeClustering(n_clusters=k, affinity=affinity, linkage=linkage).fit_predict(data)
    ari = adjusted_rand_score(y, ac)
    ami = adjusted_mutual_info_score(y, ac)
    return ari, ami


def correlation_cut(dataset, threshold):
    col_corr = set() # Set of all the names of deleted columns
    corr_matrix = dataset.corr()
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if corr_matrix.iloc[i, j] >= threshold:
                colname = corr_matrix.columns[i] # getting the name of column
                col_corr.add(colname)
                if colname in dataset.columns:
                    del dataset[colname]  # deleting the column from the dataset

    return dataset


with open('ck.pkl', 'rb') as f:
    X = pickle.load(f)

# Discovery Cohort
Xdisc = X.iloc[[i for i in range(len(X)) if 'CHD7 LOF' in X.index[i] or 'KMT2D LOF' in X.index[i] or 110 <= i < 150]]
diseases = [find_between(sample, 'blood ', ' genomic') for sample in Xdisc.index]
le = LabelEncoder()
labels = le.fit_transform(diseases)
# differentially methylated | mean: 0.9507450310749023, std: 3.3306690738754696e-16
chd7Sigs = pd.read_csv("R/chd7sigs.csv")
kmt2dSigs = pd.read_csv("R/kmt2dsigs.csv")
chd7Sigs = chd7Sigs['TargetID'].values
kmt2dSigs = kmt2dSigs['TargetID'].values
sigs = np.unique(np.append(chd7Sigs, kmt2dSigs))
XsigsDisc = filter_probes(Xdisc, sigs)


# shadi's significant probes with his betas
mbeta = pd.read_csv('m_beta_matrix.csv', header=0, index_col=0)
bestProbes = pd.read_csv('6_0005_best_probes.csv', header=0)['Probe_ID'].values
mbetaSig = np.transpose(mbeta.iloc[[i for i in range(len(mbeta)) if mbeta.index[i] in bestProbes]])

# best (ari=0.9270494745693534, ami=0.868021101679441)
hclustWrapper(mbetaSig, labels, 3, affinity='euclidean', linkage='complete')
# best (ari=0.9504501126342838, ami=0.9291765321835838)
hclustWrapper(XsigsDisc, labels, 3, affinity='euclidean', linkage='ward')

data = mbetaSig
data['cluster'] = labels
cluster = data.pop('cluster')
colDict = dict(zip(cluster.unique(), sns.hls_palette(3)))
rowCol = cluster.map(colDict)
sns.clustermap(data, method='complete', metric='euclidean', row_colors=rowCol, cmap='Blues')

data = XsigsDisc
data['cluster'] = labels
cluster = data.pop('cluster')
colDict = dict(zip(cluster.unique(), sns.hls_palette(3)))
rowCol = cluster.map(colDict)
sns.clustermap(data, method='ward', metric='euclidean', row_colors=rowCol, cmap='Blues')
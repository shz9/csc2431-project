import pandas as pd
import numpy as np
from scipy.stats import mode
import pickle
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_samples, euclidean_distances
from sklearn.preprocessing import LabelEncoder
import seaborn as sns
import matplotlib.pyplot as plt
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


def hclustWrapper(data, y, k, affinity, linkage):
    ac = AgglomerativeClustering(n_clusters=k, affinity=affinity, linkage=linkage).fit_predict(data)
    ari = adjusted_rand_score(y, ac)
    ami = adjusted_mutual_info_score(y, ac)
    return ari, ami


def kmeansWrapper(data, y, k, cores=1):
    kmeans = KMeans(n_clusters=k, n_jobs=cores).fit(data)
    ari = adjusted_rand_score(y, kmeans.labels_)
    ami = adjusted_mutual_info_score(y, kmeans.labels_)
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


def plot(df, image, dim=2):
    groups = df.groupby('disease')
    if dim == 3:
        # 3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for name, group in groups:
            ax.scatter(group.x, group.y, group.z, marker='o', label=name)

    else:
        fig, ax = plt.subplots()
        ax.margins(0.05)
        for name, group in groups:
            ax.plot(group.x, group.y, marker='o', linestyle='', ms=12, label=name)

        lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig.savefig(image, bbox_extra_artists=(lgd,), bbox_inches='tight')

    plt.show()


def kmeansPlot(data, diseases, image):
    km = KMeans(n_clusters=3, n_jobs=3).fit(data)
    chd7Class = mode(km.labels_[np.array(diseases) == 'CHD7 LOF']).mode[0]
    kmt2dClass = mode(km.labels_[np.array(diseases) == 'KMT2D LOF']).mode[0]
    distances = euclidean_distances(data, km.cluster_centers_)
    chd7Dist = max(distances[:, chd7Class]) - distances[:, chd7Class]
    kmt2dDist = max(distances[:, kmt2dClass]) - distances[:, kmt2dClass]
    distDf = pd.DataFrame(dict(x=chd7Dist, y=kmt2dDist, disease=diseases))
    plot(distDf, image)


def hclustSVM(data, n_clusters=3, n_svms=100, holdout=0.5, cores=1):
    ac = KMeans(n_clusters=n_clusters, n_init=100, max_iter=1000, n_jobs=cores).fit_predict(data)
    svc = LinearSVC(penalty='l1', loss='squared_hinge', C=1, multi_class='ovr', dual=False)

    pred = np.zeros((len(data), n_svms))
    for i in range(n_svms):
        s = random.sample(range(len(data)), int((1-holdout)*len(data)))
        svc.fit(data.iloc[s], ac[s])
        pred[:, i] = svc.predict(data)
        pred[s, i] = np.nan

    modes = mode(pred, axis=1, nan_policy='omit')[0].reshape(len(data)).astype(int)
    return modes.data


def hclustRF(data, n_clusters=3, n_rfs=100, holdout=0.5, cores=1):
    ac = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage='ward').fit_predict(data)
    rf = RandomForestClassifier(n_estimators=100, n_jobs=cores)

    pred = np.zeros((len(data), n_rfs))
    for i in range(n_rfs):
        s = random.sample(range(len(data)), int((1-holdout)*len(data)))
        rf.fit(data.iloc[s], ac[s])
        pred[:, i] = rf.predict(data)
        pred[s, i] = np.nan

    modes = mode(pred, axis=1, nan_policy='omit')[0].reshape(len(data)).astype(int)
    return modes.data


def predict_Disease(data, n_clusters=3, iters=100, stage1='kmeans', stage2='rf', stage3='kmeans', C=1, trees=100, cores=1):
    if stage1 == 'complete':
        initC = AgglomerativeClustering(n_clusters=n_clusters, linkage='complete').fit_predict(data)
    elif stage1 == 'ward':
        initC = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward').fit_predict(data)
    else:
        initC = KMeans(n_clusters=n_clusters, n_init=100, max_iter=1000, n_jobs=cores).fit_predict(data)

    bestProbes = []
    for i in range(iters):
        if stage2 == 'lr':
            model = LogisticRegression(penalty='l1', C=C)
            model.fit(data, initC)
            probes = np.unique(np.concatenate((np.where(model.coef_[0] > 0)[0],
                                               np.where(model.coef_[1] > 0)[0],
                                               np.where(model.coef_[2] > 0)[0])))
            bestProbes.extend(probes)
        elif stage2 == 'svm':
            model = LinearSVC(penalty='l1', C=C, dual=False)
            model.fit(data, initC)
            probes = np.unique(np.concatenate((np.where(model.coef_[0] > 0)[0],
                                               np.where(model.coef_[1] > 0)[0],
                                               np.where(model.coef_[2] > 0)[0])))
            bestProbes.extend(probes)
        else:
            model = RandomForestClassifier(n_estimators=trees, n_jobs=cores)
            model.fit(data, initC)
            probes = np.where(model.feature_importances_ > 0)[0]
            bestProbes.extend(probes)

    bestProbes = np.unique(bestProbes)
    bestBetas = data[data.columns[bestProbes]]

    if stage3 == 'complete':
        finalC = AgglomerativeClustering(n_clusters=n_clusters, linkage='complete').fit_predict(bestBetas)
    elif stage3 == 'ward':
        finalC = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward').fit_predict(bestBetas)
    else:
        finalC = KMeans(n_clusters=n_clusters, n_init=100, max_iter=1000, n_jobs=cores).fit_predict(bestBetas)

    return finalC, bestProbes


def hclusteRFWard(data, n_clusters=3, n_rfs=100, cores=1):
    ac = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage='ward').fit_predict(data)
    bestProbes = []
    for i in range(n_rfs):
        rf = RandomForestClassifier(n_estimators=100, n_jobs=cores)
        rf.fit(data, ac)
        probes = np.where(rf.feature_importances_ > 0)[0]
        bestProbes.extend(probes)

    bestProbes = np.unique(bestProbes)
    bestBetas = data[data.columns[bestProbes]]
    ward = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage='ward').fit_predict(bestBetas)
    return ward, bestProbes


with open('ck.pkl', 'rb') as f:
    X = pickle.load(f)
# Discovery Cohort ----------------------------------------------------------------------------------------------------
Xdisc = X.iloc[[i for i in range(len(X)) if 'CHD7 LOF' in X.index[i] or 'KMT2D LOF' in X.index[i] or 110 <= i < 150]]
diseases = [find_between(sample, 'blood ', ' genomic') for sample in Xdisc.index]
le = LabelEncoder()
labelsDisc = le.fit_transform(diseases)
# differentially methylated | mean: 0.9507450310749023, std: 3.3306690738754696e-16
chd7Sigs = pd.read_csv("R/chd7sigs.csv")
kmt2dSigs = pd.read_csv("R/kmt2dsigs.csv")
chd7Sigs = chd7Sigs['TargetID'].values
kmt2dSigs = kmt2dSigs['TargetID'].values
sigs = np.unique(np.append(chd7Sigs, kmt2dSigs))
XsigsDisc = filter_probes(Xdisc, sigs)


# shadi's significant probes with his betas
mbeta = pd.read_csv('data/discovery/m_beta_matrix.csv', header=0, index_col=0)
bestProbes = pd.read_csv('6_0005_best_probes.csv', header=0)['Probe_ID'].values
XmbetaDisc = np.transpose(mbeta.iloc[[i for i in range(len(mbeta)) if mbeta.index[i] in bestProbes]])


# ari=0.9371437539690477+/-0.04271973558404674, ami=0.8922668367721615+/-0.061524183134896764
ari = [0]*100
ami = [0]*100
for i in range(100):
    preds, topProbes = hclustRFWard(XmbetaDisc, cores=3)
    ari[i] = adjusted_rand_score(preds, labelsDisc)
    ami[i] = adjusted_mutual_info_score(preds, labelsDisc)

# Discovery + Validation Cohort ---------------------------------------------------------------------------------------
diseaseState = pd.read_csv('data/charge_kabuki/charge_kabuki_labels.csv', index_col=0, header=None)
diseaseState = diseaseState[1].values
disease2index = {'Negative': 0, 'CHARGE': 1, 'Kabuki': 2}
labelsVal = [disease2index[s] for s in diseaseState]

# shadi's significant probes with his betas
mbeta = pd.read_csv('data/validation/m_beta_matrix.csv', header=0, index_col=0)
bestProbes = pd.read_csv('6_0005_best_probes.csv', header=0)['Probe_ID'].values
XmbetaVal = np.transpose(mbeta.iloc[[i for i in range(len(mbeta)) if mbeta.index[i] in bestProbes]])

ari = [0]*100
ami = [0]*100
for i in range(100):
    preds, topProbes = hclustRFWard(XmbetaVal, cores=3)
    ari[i] = adjusted_rand_score(preds, labelsVal)
    ami[i] = adjusted_mutual_info_score(preds, labelsVal)





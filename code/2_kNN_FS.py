"""
Author: Shadi Zabad
Date: April 5th, 2018

"""

# ----------------------------
# Importing required libraries
# ----------------------------

from sklearn.neighbors import NearestNeighbors
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics import adjusted_rand_score, silhouette_score, adjusted_mutual_info_score
import pandas as pd
import numpy as np
import csv
import json
import sys
import errno
import shutil
import os
import glob
from scipy.stats import mannwhitneyu, ttest_ind
from multiprocessing import Pool, Queue

# ----------------------------
# Constants and global variables
# ----------------------------

# Constants

NUM_NN = int(sys.argv[1])
TTEST_PVAL_CUTOFF = MWUT_PVAL_CUTOFF = float(sys.argv[2])

try:
    SAVE_PROBES = bool(sys.argv[3])
except IndexError:
    SAVE_PROBES = False

if not SAVE_PROBES:
    print 'Not saving probes...'


CLUST_METHOD = 'K-means clustering'
METRIC = 'Adjusted rand'

# NUM_ITERS = 40000
# LAMBDA = 10e-4

NUM_PROCS = 100

NUM_CLUST_RUNS = 50
NUM_CLUSTERS = 3

SPECTRAL_GAMMA = .01

# Paths

RUN_PATH_SUFFIX = str(NUM_NN) + '_' + str(TTEST_PVAL_CUTOFF).replace('.', '')

BETA_MATRIX_PATH = "../output/0_processed_DNAm_data/m_beta_matrix.csv"
TARGETS_FILE_PATH = "../output/0_processed_DNAm_data/m_targets.csv"

CHD7_REF_PATH = "../metadata/Reference_CHD7_Diff_Methylated.csv"
KMT2D_REF_PATH = "../metadata/Reference_KMT2D_Diff_Methylated.csv"

U_SIG_DIR_MAIN = "../output/kNN_temp/u_tests/"
T_SIG_DIR_MAIN = "../output/kNN_temp/t_tests/"
OUTPUT_DIR_MAIN = "../output/kNN_temp/"

U_SIG_DIR = os.path.join(U_SIG_DIR_MAIN, RUN_PATH_SUFFIX + '/')
T_SIG_DIR = os.path.join(T_SIG_DIR_MAIN, RUN_PATH_SUFFIX + '/')
OUTPUT_DIR = os.path.join(OUTPUT_DIR_MAIN, RUN_PATH_SUFFIX + '/')


# ----------------------------
# Auxiliary functions
# ----------------------------


def write_feature_list_to_csv(file_name, feat_data):

    with open(file_name, "wb") as outf:

        csv_writer = csv.writer(outf)
        csv_writer.writerow(['ProbeID', 'pval'])

        for sfeat in feat_data:
            csv_writer.writerow(sfeat)


def read_and_concat_files(dir_path):

    df_list = []

    for fname in glob.glob(dir_path + "/*.csv"):
        df_list.append(pd.read_csv(fname, index_col=None, header=0))

    final_df = pd.concat(df_list)
    final_df.columns = ['Probe_ID', 'pval']

    return final_df


def bootstrap_clustering(filt_data, method='kmeans'):

    ar_scores = []
    ami_scores = []
    silh_scores = []

    for ridx in range(NUM_CLUST_RUNS):

	if method == 'kmeans':
        clust_data = KMeans(n_clusters=NUM_CLUSTERS, n_jobs=-1).fit(filt_data)
    elif method == 'spectral':
        clust_data = SpectralClustering(n_clusters=NUM_CLUSTERS, gamma=SPECTRAL_GAMMA, n_jobs=-1).fit(filt_data)
	else:
        raise Exception("Clustering method not supported.")

    ar_scores.append(adjusted_rand_score(exp_targets['Disease_State'], clust_data.labels_))
    ami_scores.append(adjusted_mutual_info_score(exp_targets['Disease_State'], clust_data.labels_))
    silh_scores.append(silhouette_score(filt_data, clust_data.labels_))

    return {'Adjusted rand': {'Mean': np.mean(ar_scores), 'std': np.std(ar_scores)},
            'Mutual information': {'Mean': np.mean(ami_scores), 'std': np.std(ami_scores)},
            'Silhouette': {'Mean': np.mean(silh_scores), 'std': np.std(silh_scores)}}



def eval_probe_selection(probe_sel_criteria, probe_ids):

    chd7_overlap = len(chd7_ref[chd7_ref['DIFF_METH_PROBES'].isin(probe_ids['Probe_ID'])])
    kmt2d_overlap = len(kmt2d_ref[kmt2d_ref['DIFF_METH_PROBES'].isin(probe_ids['Probe_ID'])])

    filtered_data = beta_mat[list(probe_ids['Probe_ID'])]

    score_dict = {
        'Probe selection criteria': probe_sel_criteria,
        'Number of selected probes': len(probe_ids['Probe_ID']),
        'Overlap with reference probes': {
            'CHD7': chd7_overlap,
            'KMT2D': kmt2d_overlap
        },
        'K-means clustering': bootstrap_clustering(filtered_data),
        'Spectral clustering': bootstrap_clustering(filtered_data, method='spectral')
    }

    return score_dict


# ----------------------------
# Workflow steps
# ----------------------------

# In what follows, I'll implement (with some modifications) the algorithm described in:
#
# Rethinking Unsupervised Feature Selection: From Pseudo Labels to Pseudo Must-links
# Wei et al. 2017
#
# https://www.cs.uic.edu/~bcao1/doc/ecmlpkdd17b.pdf


# Create directories if they don't exist

for cdir in (OUTPUT_DIR, U_SIG_DIR, T_SIG_DIR):
    try:
	os.makedirs(OUTPUT_DIR)
    except OSError as e:
	if e.errno != errno.EEXIST:
        sys.exit(str(e))

# Remove already existing files:

try:
    shutil.rmtree(U_SIG_DIR, ignore_errors=True)
    os.makedirs(U_SIG_DIR)

    shutil.rmtree(T_SIG_DIR, ignore_errors=True)
    os.makedirs(T_SIG_DIR)

except Exception, e:
    if e.errno != errno.EEXIST:
        sys.exit(str(e))


# Read the targets file:

exp_targets = pd.read_csv(TARGETS_FILE_PATH, header=0, index_col=0)

# Read the beta matrix:

beta_mat = pd.read_csv(BETA_MATRIX_PATH, header=0, index_col=0)
beta_mat = beta_mat.T  # convert to matrix of dims (n_patients, n_probes)

probe_ids = beta_mat.columns

# Find nearest neighbors
nbrs = NearestNeighbors(n_neighbors=NUM_NN, n_jobs=-1).fit(beta_mat)

# Generate the neighbor graph
neighbor_graph = nbrs.kneighbors_graph(beta_mat).toarray()

print "Neighbor graph dimensions: ", neighbor_graph.shape

# Create matrices for similar and dissimilar pairs.
# These matrices will be of dimensions (n_pairs, n_probes)
# To combine the probe values for a pair, I'll follow Wei et al. and multiply the values.


must_link = np.empty((0, len(probe_ids)))
cannot_link = np.empty((0, len(probe_ids)))

similar_pairs = []
dissimilar_pairs = []

for subj1 in range(neighbor_graph.shape[0]):
    for subj2 in range(subj1 + 1, neighbor_graph.shape[0]):
        if neighbor_graph[subj1, subj2] == 1 or neighbor_graph[subj2, subj1] == 1:
            must_link = np.append(must_link,
                                  np.array([np.multiply(beta_mat.iloc[subj1, ], beta_mat.iloc[subj2, ])]),
                                  axis=0)
            similar_pairs.append((subj1, subj2))
        else:
            cannot_link = np.append(cannot_link,
                                    np.array([np.multiply(beta_mat.iloc[subj1, ], beta_mat.iloc[subj2, ])]),
                                    axis=0)
            dissimilar_pairs.append((subj1, subj2))

# ----------------------------
# * * * Hypothesis testing frameworks * * *
# ----------------------------

# Here, I'll experiment with 2 hypothesis testing frameworks to see if
# the values of the features in similar pairs are distinguishable from those
# in dissimilar pairs.

# Apply the Mann-Whitney U and Student-t tests:


def apply_tests(pidx):
    # MWUT:
    u, u_pval = mannwhitneyu(must_link[:, pidx], cannot_link[:, pidx])

    if u_pval < MWUT_PVAL_CUTOFF:
        write_feature_list_to_csv(U_SIG_DIR + str(pidx) + "_u_sig_features.csv", [(probe_ids[pidx], u_pval)])

    # Student-t:
    t, t_pval = ttest_ind(must_link[:, pidx], cannot_link[:, pidx])

    if t_pval < TTEST_PVAL_CUTOFF:
        write_feature_list_to_csv(T_SIG_DIR + str(pidx) + "_t_sig_features.csv", [(probe_ids[pidx], t_pval)])


def launch_pool_process():
    t_pool = Pool(NUM_PROCS)
    t_pool.map(apply_tests, list(range(len(probe_ids))))

    t_pool.close()
    t_pool.join()

if __name__ == '__main__':
    print "Launching pool..."
    launch_pool_process()


u_sig_probes = read_and_concat_files(U_SIG_DIR)
t_sig_probes = read_and_concat_files(T_SIG_DIR)

intersect_probes = u_sig_probes[u_sig_probes['Probe_ID'].isin(t_sig_probes['Probe_ID'])]

u_sig_probes.to_csv(OUTPUT_DIR + "u_test_results.csv")
t_sig_probes.to_csv(OUTPUT_DIR + "t_test_results.csv")

chd7_ref = pd.read_csv(CHD7_REF_PATH)
kmt2d_ref = pd.read_csv(KMT2D_REF_PATH)

num_chd7_ref = len(chd7_ref)
num_kmt2d_ref = len(kmt2d_ref)

result_set = {
    'Nearest neighbors': NUM_NN,
    'U-test pval': MWUT_PVAL_CUTOFF,
    't-test pval': TTEST_PVAL_CUTOFF,
    'Number of reference probes': {
	'CHD7': num_chd7_ref,
	'KMT2D': num_kmt2d_ref
    },
    'Evaluation': []
}

probe_sel_dict = dict(zip(['t-test probes', 'u-test probes', 'intersect probes'],
                          [t_sig_probes, u_sig_probes, intersect_probes]))

for s_name, pids in probe_sel_dict.iteritems():
    result_set['Evaluation'].append(eval_probe_selection(s_name, pids))

result_set['Evaluation'] = sorted(result_set['Evaluation'],
                                  key=lambda x:x[CLUST_METHOD][METRIC]['Mean'],
                                  reverse=True)

with open(OUTPUT_DIR + 'result.json', 'wb') as outf:
    json.dump(result_set, outf)

if SAVE_PROBES:
    best_performing_probes = probe_sel_dict[result_set['Evaluation'][0]['Probe selection criteria']]
    best_performing_probes.to_csv(OUTPUT_DIR + 'best_performing_probes.csv')

# ----------------------------
# * * * Classification frameworks * * *
# Code below works, but it's not giving good results.. ignoring for now.
# ----------------------------

"""
def subgradient_descent_clf():

    w_vec = np.zeros(len(probe_ids))
    w_mat = np.zeros((len(probe_ids), len(probe_ids)))

    for i in range(NUM_ITERS):

        if np.random.rand() > .5:
            s_pair = similar_pairs[np.random.choice(len(similar_pairs))]
            l_ij = 1
        else:
            s_pair = dissimilar_pairs[np.random.choice(len(dissimilar_pairs))]
            l_ij = -1

	np.fill_diagonal(w_mat, w_vec)
	s_ij = np.dot(np.dot(beta_mat.iloc[s_pair[0], ].T, w_mat), beta_mat.iloc[s_pair[1], ])

	w_vec += LAMBDA * np.sign(w_vec)

	if s_ij * l_ij < 1:
        w_vec += l_ij * np.multiply(beta_mat.iloc[s_pair[0], ], beta_mat.iloc[s_pair[1], ])


    return w_vec


clf_weights = subgradient_descent_clf()
probe_weights = zip(probe_ids, clf_weights)
probe_weights = sorted(probe_weights, key=lambda x: np.abs(x[1]), reverse=True)

print probe_weights[:100]

sorted_probes = [p[0] for p in probe_weights]

print "For CHD7 (" + str(num_chd7_ref) + "):"
print "---------"
print "Overlap in top 500:", len(chd7_ref[chd7_ref['DIFF_METH_PROBES'].isin(sorted_probes[:500])])
print "Overlap in top 1000:", len(chd7_ref[chd7_ref['DIFF_METH_PROBES'].isin(sorted_probes[:1000])])
print "Overlap in top 2000:", len(chd7_ref[chd7_ref['DIFF_METH_PROBES'].isin(sorted_probes[:2000])])
print "Overlap in top 5000:", len(chd7_ref[chd7_ref['DIFF_METH_PROBES'].isin(sorted_probes[:5000])])

print "For KMT2D (" + str(num_kmt2d_ref) + "):"
print "---------"
print "Overlap in top 500:", len(kmt2d_ref[kmt2d_ref['DIFF_METH_PROBES'].isin(sorted_probes[:500])])
print "Overlap in top 1000:", len(kmt2d_ref[kmt2d_ref['DIFF_METH_PROBES'].isin(sorted_probes[:1000])])
print "Overlap in top 2000:", len(kmt2d_ref[kmt2d_ref['DIFF_METH_PROBES'].isin(sorted_probes[:2000])])
print "Overlap in top 5000:", len(kmt2d_ref[kmt2d_ref['DIFF_METH_PROBES'].isin(sorted_probes[:5000])])
"""


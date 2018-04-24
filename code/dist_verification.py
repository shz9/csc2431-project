"""
Author: Shadi Zabad
Date: April 2018

This script is used to test the robustness of various
distance metrics in separating patients based on
their DNA methylation profiles (Beta values).
"""

# ----------------------------
# Importing required libraries
# ----------------------------

import pandas as pd
import numpy as np

from sklearn.neighbors import DistanceMetric

import matplotlib
matplotlib.use('Agg')


# ----------------------------
# Constants and global variables
# ----------------------------

BETA_MATRIX_PATH = "../output/0_processed_DNAm_data/m_beta_matrix.csv"
TARGETS_FILE_PATH = "../output/0_processed_DNAm_data/m_targets.csv"

# ----------------------------
# Auxiliary functions
# ----------------------------

def plot_distance_matrix(metric, dist_mat):

    plt.matshow(dist_mat)
    plt.savefig("../plots/" + metric + "_distance_mat.png")
    plt.close()


def plot_sorted_label_matrix(metric, dist_mat, true_labels, disease_labels=True):

    li = []

    if disease_labels:
        lab = {
            'CHARGE': 0.0,
            'Kabuki': 1.0,
            'Control': 2.0
        }
    else:
        lab = {
            'CHD7_LOF_discovery_cohort': 0.0,
            'KMT2D_LOF_discovery_cohort': 1.0,
            'Control_for_CHD7_LOF_discovery_cohort': 2.0,
            'Control_for_KMT2D_LOF_discovery_cohort': 3.0
        }

    for subj in range(len(dist_mat)):
        labels = [lab[true_labels[lb]] for lb in np.argsort(dist_mat[subj])]
        li.append(labels)

    arr = np.array(li)

    if disease_labels:
        suff = 'label'
    else:
        suff = 'cohort'

    plt.matshow(arr)
    plt.xticks(np.arange(-.5, len(arr)), [])
    plt.yticks(np.arange(-.5, len(arr)), [])
    plt.grid(which='major', color='w', linestyle='-', linewidth=.5)
    plt.savefig("../plots/" + metric + "_" + suff + "_mat.png")
    plt.close()


# ----------------------------
# Workflow steps
# ----------------------------

# Read the targets file:

exp_targets = pd.read_csv(TARGETS_FILE_PATH, header=0, index_col=0)

# Read the beta matrix:

beta_mat = pd.read_csv(BETA_MATRIX_PATH, header=0, index_col=0)
beta_mat = beta_mat.T  # convert to matrix of dims (n_patients, n_probes)

distances_to_test = ['euclidean', 'manhattan', 'chebyshev',
                     'minkowski'] # , 'mahalanobis']

dist_dict = dict()

for dist_name in distances_to_test:

    if dist_name == 'mahalanobis':
        dist = DistanceMetric.get_metric(dist_name, V=np.cov(beta_mat,
                                                             rowvar=False))
    else:
        dist = DistanceMetric.get_metric(dist_name)

    dist_dict[dist_name] = dist.pairwise(beta_mat)


true_labels = exp_targets['Disease_State']
NN = 6

for metric in dist_dict.keys():
    plot_distance_matrix(metric, dist_dict[metric])
    plot_sorted_label_matrix(metric, dist_dict[metric], true_labels)
    plot_sorted_label_matrix(metric, dist_dict[metric],
                             exp_targets['Sample_Group'],
                             disease_labels=False)
    print '>>>', metric
    for subj in range(len(true_labels)):
        subj_label = true_labels[subj]
        closest = np.argsort(dist_dict[metric][subj])[: NN + 1]
        closest_labels = [true_labels[i] for i in closest if i != subj]
        print exp_targets.index[subj], float(closest_labels.count(subj_label)) / NN


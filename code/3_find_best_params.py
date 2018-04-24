"""
Author: Shadi Zabad
Date: April 18th, 2018
"""

# ----------------------------
# Importing required libraries
# ----------------------------

import sys
import os
import fnmatch
import json
import pprint

# ----------------------------
# Constants and global variables
# ----------------------------

RESULT_SET_DIR = "../output/kNN_temp/"
OUTPUT_FILE_PATTERN = "result.json"

METRICS = {
    'AR': 'Adjusted rand',
    'MI': 'Mutual information',
    'SI': 'Silhouette'
}

CLUSTERING_METHODS = {
    'KM': 'K-means clustering',
    'SP': 'Spectral clustering'
}

# Select the quality metric to sort results by:
try:
    METRIC = sys.argv[1]
except IndexError:
    METRIC = 'AR'


if METRIC not in METRICS.keys():
    raise Exception('Metric {} is not supported.'.format(METRIC))
else:
    METRIC = METRICS[METRIC]


# Select the clustering method to sort results by:
try:
    CLUST_METHOD = sys.argv[2]
except IndexError:
    CLUST_METHOD = 'KM'


if CLUST_METHOD not in CLUSTERING_METHODS.keys():
    raise Exception('Clustering method {} is not supported.'.format(CLUST_METHOD))
else:
    CLUST_METHOD = CLUSTERING_METHODS[CLUST_METHOD]

# ----------------------------
# Workflow steps
# ----------------------------

res_list = []

for root, dirnames, filenames in os.walk(RESULT_SET_DIR):
    for filename in fnmatch.filter(filenames, OUTPUT_FILE_PATTERN):
	with open(os.path.join(root, filename), 'rb') as rf:
	    res_list.append(json.load(rf))

for idx in range(len(res_list)):
    res_list[idx]['Evaluation'] = sorted(res_list[idx]['Evaluation'], 
					 key=lambda x:x[CLUST_METHOD][METRIC]['Mean'],
					 reverse=True)

res_list = sorted(res_list, 
		  key=lambda x:x['Evaluation'][0][CLUST_METHOD][METRIC]['Mean'], 
		  reverse=True)

print '##############################################'
print 'Best NN params, in order:', [(p['Nearest neighbors'], p['U-test pval']) for p in res_list]
print '##############################################'


for elem in res_list:
    pprint.pprint(elem)
    print '-------------------------------'



# Plot neighborhood diversity index for UMAP
# Adapted from Arjun Rao's code: 
#    https://github.com/UCSF-DSCOLAB/combes_et_al_COVID_2020/blob/master/auxiliary_scripts/batch_effect_plot.py 
# See Extended Data Fig 1b. (https://www.nature.com/articles/s41586-021-03234-7/figures/5) 
#   from Combes et al Nature 2021 for an example
#
# The neighborhood diversity index refers to the number of a cell's neighbors from other groups.
# This is particularly useful for visualizing datasets where there are a lot batches or groups.
#
# To use the script, 
# (1) first generate a tsv with the following columns: idx (rownumber), LIBRARY (aka batch), UMAP_1, UMAP_2
#    For large datasets, it is recommended to downsample the number of cells to 100-200k 
#    to decrease computation time.
# (2) then run:
#       python3 plot_div_index.py obj_dir obj_name
#  obj_dir is the directory the file is located in
#  obj_name is the name of the .tsv file minus its suffix
#
# This script will then output the batch plots as "<obj_name>_batch_plot.pdf"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from collections import Counter
from sklearn.neighbors import NearestNeighbors

import os
import sys


obj_dir = sys.argv[1]
obj_name = sys.argv[2]
os.chdir(obj_dir)

x = pd.read_csv('{}.tsv'.format(obj_name), sep='\t', header=0, index_col=0)
x = x[['LIBRARY', 'UMAP_1', 'UMAP_2']].copy()

y = np.array(list(zip(x['UMAP_1'], x['UMAP_2'])))
nbrs = NearestNeighbors(n_neighbors=int(x.shape[0]**0.5), algorithm='ball_tree').fit(y)
distances, indices = nbrs.kneighbors(y)

threshold = float(np.quantile(distances, q=[0.9]))
libraries = np.array(list(x['LIBRARY']))
library_avgs = Counter(x['LIBRARY'])
library_avgs = {k: v/x.shape[0] for k, v in library_avgs.items()}

def dscore(nbh, avgs):
  temp = Counter(nbh)
  temp = Counter({k: v/len(nbh) for k, v in temp.items()})
  return 0.5 * sum([abs(avgs[k] - temp[k]) for k in avgs])

x['ND'] = [len(set(libraries[i[d<=threshold]])) for d, i in zip(distances, indices)]
plt.figure(figsize=(10, 10))
x.plot.scatter(x='UMAP_1', y='UMAP_2', c='ND', s=0.1, vmin=0, vmax=len(set(libraries)), cmap="mako") 
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.axis('off')
plt.savefig('{}_batch_plot.pdf'.format(obj_name), bbox_inches='tight', pad_inches=0.0)
plt.close()

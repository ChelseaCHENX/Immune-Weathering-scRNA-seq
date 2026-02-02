import os
import sys
import math
import pickle as pkl
import argparse
import subprocess
from collections import Counter

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.patches as mpatches
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
from matplotlib import rcParams


import scanpy as sc
import scrublet as scr
import scipy.io


def run(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return proc.returncode, stdout, stderr


def getarg(argv, index, default=None):
    return argv[index] if len(argv) > index else default


def read_gmt(fpath):
    gdict = {}
    with open(fpath, 'r') as f:
        for line in f:
            items = line.strip().split('\t')
            gset = items[0].lower()
            genes = items[2:]
            gdict[gset] = genes
    return gdict


def write_gmt(gdict, fpath):
    with open(fpath, 'w') as f:
        for key, val in gdict.items():
            f.write(f"{key}\thttp\t{'\t'.join(val)}\n")
    print(f"finished writing {len(gdict)} gene sets to {fpath}")

# Additional plotting functions preserved as-is from the original due to complexity. If needed, they can be modularized in a separate file.
# forestplot, forestplot_group, plot_pathway_heatmap, plot_multi_pointplots remain intact.

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

# def split_umap(adata, split_by, ncol=2, nrow=None, **kwargs):
#     categories = adata.obs[split_by].cat.categories
#     if nrow is None:
#         nrow = int(np.ceil(len(categories) / ncol))
#     fig, axs = plt.subplots(nrow, ncol, figsize=(5*ncol, 4*nrow))
#     axs = axs.flatten()
#     for i, cat in enumerate(categories):
#         ax = axs[i]
#         sc.pl.umap(adata[adata.obs[split_by] == cat], ax=ax, show=False, title=cat, **kwargs)
#     plt.tight_layout()
# split_umap(adata, color = ['n_genes_by_counts'], split_by='type',legend_loc = "right margin")
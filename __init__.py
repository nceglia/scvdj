import scanpy as sc
import numpy
import sklearn
import pandas
import textdistance
import collections
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
from matplotlib.colors import ListedColormaps
import operator
import math

def read_10x_vdj(adata, directory):
    return adata

def tcr_plot(adata, reduction="umap", top_n=10, filename="tcr_plot.png"):
    df = adata.obs
    dft = df[df["cdr3s_aa"].notnull()]
    plt.figure(figsize = (7, 7))
    clonotype_counts = collections.defaultdict(int)
    for clonotype in df['cdr3s_aa']:
        clonotype_counts[clonotype] += 1
    top_clonotypes = sorted(clonotype_counts.items(), key=operator.itemgetter(1),reverse=True)
    top_clonotypes = [x[0] for x in top_clonotypes[:top_n]]
    ax1 = plt.subplot(1,1,1)
    sizes = []
    clonotype_labels = []
    cluster_labels = []
    for clonotype in df['cdr3s_aa']:
        clonotype = str(clonotype)
        if clonotype not in top_clonotypes or clonotype == "Other" or clonotype == "nan":
            sizes.append(1) 
            clonotype_labels.append("_Other") 
        else:
            sizes.append(12) 
            clonotype_labels.append(str(clonotype) + " {}".format(clonotype_counts[clonotype]))
    df["Size"] = sizes
    df["TCR Sequence"] = clonotype_labels
    colors = ["#75715E", "#fd971e","#66D9EF","#E6DB74","#7FFFD4","#E69F66","#AE81FF","#000000","#f92572","#A6E22E","#f47d7c","#0a012d","#3aaea3","#cf4301","#999999"]
    customPalette = sns.set_palette(sns.color_palette(colors))
    sns.scatterplot(data=df,x="UMAP_1", y="UMAP_2", hue='TCR Sequence', ax=ax1, alpha=0.7,linewidth=0.00,size=sizes,palette=customPalette)
    ax1.set_xlabel('UMAP-1')
    ax1.set_ylabel('UMAP-2')
    ax1.xaxis.set_ticklabels([])
    ax1.yaxis.set_ticklabels([])
    ax1.xaxis.set_ticks([])
    ax1.yaxis.set_ticks([])
    ax1.set_title("Top TCR Sequence")
    h,l = ax1.get_legend_handles_labels()
    ax1.legend(h[:top_n], l[:top_n], borderaxespad=0.,fontsize='9')
    plt.tight_layout()
    plt.savefig(filename,dpi=200)

def distance_cluster(adata, level="nt", chain="trb", edit_distance=1):
    pass

def nt_convergence(adata):
    pass
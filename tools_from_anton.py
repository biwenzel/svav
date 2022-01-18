import sys
from matplotlib import pyplot as plt
from Bio import Align
from Bio.Emboss.Applications import NeedleCommandline
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
import numpy as np
from time import time


def prepare_dist_matrix(seqs):
    n = len(seqs)
    dist_mat = np.zeros((n, n))

    aligner = Align.PairwiseAligner(mode="global", match_score=2, open_gap_score=-0.5,extend_gap_score=-0.1, mismatch_score=-0.1)

    for i in range(n):
        for j in range(i, n):
            if seqs[i]=='' and seqs[j]=='':
                result = 0
            elif seqs[j]=='':
                result = (-0.8+(len(seqs[i])-1)*-0.2)/len(seqs[i])
            elif seqs[i]=='':
                result = (-0.8+(len(seqs[j])-1)*-0.2)/len(seqs[j])
            else:
                alignment = aligner.align(seqs[i], seqs[j])[0]
                al_len = len(format(alignment).split("\n")[0])
                result = alignment.score / al_len
            dist_mat[j, i] = dist_mat[i, j] = result

    #print(dist_mat.tolist())
    return dist_mat


def pca_on_read_distance(seqs):
    starttime = time()
    dist_mat = prepare_dist_matrix(seqs)
    print(f"Distances for PCA calculated (in {time()-starttime:.2f}s)...")

    pca = PCA(n_components=2)
    res = pca.fit_transform(dist_mat)
    return res

def cluster_reads(data, eps=0.15, min_samples=5, normalize_scales=True):
    data_norm = data.copy()
    if normalize_scales:
        for axis in range(data.shape[1]):
            data_norm[:, axis] /= (data_norm[:, axis].max() - data_norm[:, axis].min())
    res = DBSCAN(eps=eps, min_samples=min_samples).fit(data_norm)
    return res

def plot_pca_result(data, labels=None):
    """ labels are going to be interpreted the following way:
        the first n - (len(labels)-1) points are labeled with the first label
        then, all the remaining dots get their own label """
    x, y = zip(data.T)
    x, y = x[0], y[0]
    if labels is None or len(labels) == 0:
        plt.scatter(x, y)
    else:
        n = len(x)
        read_count = n - (len(labels)-1)
        plt.scatter(x[:read_count], y[:read_count])
        for i in range(read_count, n):
            plt.scatter(x[i], y[i])
        plt.legend(labels)

def plot_with_clustering(data, clustering, additional_labels=None):
    x, y = zip(data.T)
    x, y = x[0], y[0]

    names = []
    cluster_colors = plt.cm.Set3(np.linspace(0, 1, int(clustering.labels_.max())+1))
    for i, c in enumerate(cluster_colors):
        idcs = clustering.labels_ == i
        plt.scatter(x[idcs], y[idcs], color=c)
        names.append(f"Read Cluster {i}")

    if -1. in clustering.labels_:
        idcs = clustering.labels_ == -1.
        plt.scatter(x[idcs], y[idcs], color="black")
        names.append("Read outliers")

    if not additional_labels is None:
        n_add = len(additional_labels)
        for i, _ in enumerate(additional_labels):
            idx = -n_add+i
            plt.scatter(x[idx], y[idx])

        names += additional_labels

    plt.legend(names)

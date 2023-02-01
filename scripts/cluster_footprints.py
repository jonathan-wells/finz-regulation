#!/usr/bin/env python3

import re
import numpy as np
import markov_clustering as mc
import networkx as nx
from Bio import SeqIO

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return set(kmers)

def get_footprint(znfrecord):
    # znfdoms =  re.findall(r'C..C(.{12})H...H', str(znfrecord.seq))
    # return set(''.join((i[-1], i[-4], i[-5], i[-7])) for i in znfdoms)
    return build_kmers(znfrecord.seq, 10)

def load_znfs(filename):
    return list(SeqIO.parse(filename, 'fasta'))

def jaccard(set1, set2):
    return len(set1.intersection(set2))/len(set1.union(set2))

def generate_adjacency_graph(znfs):
    footprints = [(znfrecord.name, get_footprint(znfrecord)) for znfrecord in znfs]
    adjacency_list = []
    for i, (znf1, fp1) in enumerate(footprints):
        for j, (znf2, fp2) in enumerate(footprints):
            if j >= i:
                break
            adjacency_list.append((znf1, znf2, jaccard(fp1, fp2)))
    return adjacency_list

def cluster_matrix(adjacency_list):
    with open('../data/finz_clusters.txt') as infile:
        finz_type = {i.split()[0]: i.split()[1] for i in infile}
    footprint_graph = nx.Graph()
    footprint_graph.add_weighted_edges_from(adjacency_list)
    labels = [name[:-3] for name in footprint_graph.nodes]
    matrix = nx.to_scipy_sparse_array(footprint_graph)
    result = mc.run_mcl(matrix, inflation=4.0)
    clusters = []

    for cluster in mc.get_clusters(result):
        c = [finz_type[labels[i]] for i in list(cluster)]
        c = [i for i in c if i != 'NA']
        clusters.append(c)
        if len(c) != 0:
            print(clusters[-1])
    print(len(clusters))

if __name__ == "__main__":
    znfs = load_znfs('../data/seqs/Danio_rerio_hiqual_finz.aa.fa')
    adjacency_list = generate_adjacency_graph(znfs)
    cluster_matrix(adjacency_list)

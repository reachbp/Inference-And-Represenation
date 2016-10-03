
from __future__ import division
__author__ = 'bharathipriyaa'

import math
import scipy
import numpy as np
import pickle
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
input = "../data/chowliu-input.txt"
namesFile = "../data/names.txt"

"""
Load marginals methods calculated the marginal probabilities for
all x_i = count(x_i)/N
"""
def load_marginals():
    print "Loading the input file"
    df = np.loadtxt(input, dtype=int)
    edgeMatrix = df.sum(axis=0)
    N,M =  df.shape[0],df.shape[1]
    print "Calculated marginal probabilities"
    mprob = [x / N for x in edgeMatrix]
    pickle.dump(mprob, open( "marg.p", "wb" ) )

"""
Load data does the following :
1) Calculate the adjacency matrix for each node x_i, x_j.
2) Store the node potentials (in this case count(x_i,x_j) for each possible x_i,x_j in a 3d array
3) The edge weights are calculated as per Chow-liu algorithms
"""
def load_data():
    print "Loading the input file"
    df = np.loadtxt(input, dtype=int)
    edgeMatrix = df.sum(axis=0)
    N,M =  df.shape[0],df.shape[1]
    # Load marginal probabilities
    mprob = pickle.load(open("marg.p", "rb"))
    adjMatrix = np.zeros((M,M,4))
    print "Calculating Adjacency matrix"
    for i in range(N - 1):
        print i
        for j in range(M - 1):
            # f counters for (00) (01) (10) (11)
            for k in range(j + 1, M - 1):
                index1 = df[i][j] * 2 + df[i][k]
                adjMatrix[j][k][index1] += 1
                index2 = df[i][j]  + df[i][k] * 2
                adjMatrix[k][j][index2] += 1
    pickle.dump(adjMatrix, open( "adj1.p", "wb" ) )

    adjMatrix = pickle.load(open("adj1.p", "rb"))
    print "Loaded the adjMatrix"
    edgeWeight = np.zeros((M,M))
    print "Calculating edge weights"
    for i in range(M - 1):
        for j in range(i + 1, M - 1):
            for index in range(3):
                    p_xi = mprob[i] if index / 2 == 1 else 1 - mprob[i]
                    p_xj = mprob[j] if index % 2 == 1 else 1 - mprob[j]
                    p_xixj = adjMatrix[i][j][index]/ N
                    if(p_xixj != 0) :
                        edgeWeight[i][j] += (-1) * p_xixj * (math.log(p_xixj) - math.log(p_xi*p_xj))
            edgeWeight[j][i] = edgeWeight[i][j]
    print "Size of adjacency matrix {}".format(adjMatrix.shape)
    pickle.dump(edgeWeight, open( "edge.p", "wb" ) )

"""
With edge weights inverted, Run the minimum spanning tree algo.
This gives the maximum spanning tree.
Output the resulting edges of MST to .dot file

"""
def run_mst() :
    edgeWeight = pickle.load(open( "edge.p", "rb"))
    print " Determine Maximum spanning tree"
    X = csr_matrix(edgeWeight)
    names = np.loadtxt(namesFile, dtype = basestring)
    Tcsr = minimum_spanning_tree(X)
    cx = scipy.sparse.coo_matrix(Tcsr)
    print cx.shape
    print "Printing maximum spanning tree"
    str ="graph {"
    print len(cx.row)
    for i,j,v in zip(cx.row, cx.col, cx.data):
        str += names[i] + "--" + names[j] + ";"
        print "(%s, %s), %f" % (names[i],names[j],-v)
    str += "}"
    f = open('../data/graph.dot', 'w')
    f.write(str)
    pickle.dump(cx, open("cx.p", "wb"))
    print "Dumped maximum spanning tree"

"""
Print the resulting Markov network as a UAI format
"""
def print_uai():
    # Load mst_matrix
    cx = pickle.load(open("cx.p", "rb"))
    # Load marginal probabilities
    mprob = pickle.load(open("marg.p", "rb"))
    # Load adjacency matrix
    adjMatrix = pickle.load(open("adj1.p", "rb"))
    print "Pring UAI format"
    str = "MARKOV" + "\n"
    str += "{}\n".format(cx.shape[0])
    for i in range(cx.shape[0]):
        str += "{} ".format(2)
    str += "\n"
    # Add number of cliques to UAI
    for i,j,v in zip(cx.row, cx.col, cx.data):
        str += "2 {} {} \n".format(i, j)
    # Add node potentials
    for i in range(cx.shape[0]):
        str += "2 \n"
        str += "{} {}\n".format(1-mprob[i], mprob[i])
        str += "\n"
    # Add edge potentials
    for i,j,v in zip(cx.row, cx.col, cx.data):
        total = 0;
        for index in range(3):
            total += adjMatrix[i][j][index]
        str += "4 \n"
        for index_1 in {0,1}:
            for indeX_2 in {0,1}:
                str += "{} ".format(adjMatrix[i][j][index_1*2 + indeX_2]/total)
            str += "\n"
        str += "\n"
    f = open('../data/uai.dot', 'w')
    f.write(str)

if __name__ == '__main__':
    #load_marginals()
    #load_names()
    #load_data()
    #calculate_edge()
    #run_mst()
    print_uai()
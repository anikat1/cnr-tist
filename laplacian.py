# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import networkx as nx
import pandas as pd
import sys
from scipy.io import savemat

if __name__ == '__main__':
	fileName = sys.argv[1]#'data/Harvey_60min_sample.csv';
	adjlist = sys.argv[2]#'data/Harvey.adjlist';
	savefile= sys.argv[3] #'data/Harvey_laplace.mat'
	data = pd.read_csv(fileName)
	counties = list(data.columns.values)
	#print(len(counties))
	G = nx.read_adjlist(adjlist,delimiter='\t')
	L = nx.laplacian_matrix(G,counties)
	print(L)
	savemat(savefile,{'L':L})


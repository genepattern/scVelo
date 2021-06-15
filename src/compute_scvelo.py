import os, sys
import subprocess

# subprocess.check_call(['apt-get', 'update'])
# subprocess.check_call(['apt-get', 'install', '-y', 'python3-pip'])
# 
# import pkg_resources
# 
# required = {'llvmlite','numpy','anndata','scipy','pandas','scanpy','python-igraph','matplotlib', 'scvelo', 'louvain', 'pybind11', 'hnswlib'}
# installed = {pkg.key for pkg in pkg_resources.working_set}
# missing = required - installed
# 
# if missing:
#     # implement pip as a subprocess:
#     subprocess.check_call([sys.executable, '-m', 'pip', 'install',*missing])

from optparse import OptionParser
import argparse
import shutil

import os,sys
import anndata
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
import matplotlib
import igraph


def main():
	usage="%prog [options]" + "\n"
	ap = argparse.ArgumentParser()
	ap.add_argument("-i","--input-file",action="store",dest="input_file",help="h5ad (anndata) file.")
	ap.add_argument("-m","--markers",action="store",dest="markers",help="A list of marker genes")
	ap.add_argument("-s","--shared",action="store",dest="minshared",help="Filter genes by minimum shared counts")
	ap.add_argument("-t","--top",action="store",dest="topgenes",help="Top Genes for Velocity Computation")
	ap.add_argument("-e","--embedding",action="store",dest="embedding",help="Dataset was processed with umap or tsne embedding")
	ap.add_argument("-o","--out",action="store",dest="output",help="Output file basename")
	ap.add_argument("-j","--cpu",action="store",dest="ncores",help="CPU cores to use for transition dynamics calculation")

	options = ap.parse_args()

	adata=anndata.read_h5ad(options.input_file)
	
	if bool(options.markers):
		with open(options.markers) as f:
			markergenes = f.read().splitlines()
		markergenes = list(set([sub.replace('-I', '') for sub in markergenes]))


	scv.pp.filter_and_normalize(adata, min_shared_counts=int(options.minshared), n_top_genes=int(options.topgenes), enforce=True)
	scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
	scv.tl.recover_dynamics(adata)
	scv.tl.velocity(adata, mode = 'dynamical')
	scv.tl.velocity_graph(adata)
	scv.tl.umap(adata)
	scv.tl.louvain(adata)
	scv.pl.velocity_embedding_stream(adata, basis=options.embedding,save="embedding")
	ad.AnnData.write(adata, options.output + "_graph_result.h5ad")

	if bool(options.markers):
		for i in markergenes:
			sc.pl.umap(adata, color=[markergenes[i]], save="_"+options.output+"_"+markergenes[i])

if __name__ == '__main__':
	main()

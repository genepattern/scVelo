import os, sys
import subprocess
# 
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
	ap.add_argument("-v","--hvg",action="store",dest="hvg",help="Compute highly_variable_genes")
	ap.add_argument("-e","--embedding",action="store",dest="embedding",help="Dataset was processed with umap or tsne embedding")
	ap.add_argument("-o","--out",action="store",dest="output",help="Output file basename")
	ap.add_argument("-j","--cpu",action="store",dest="ncores",help="CPU cores to use for transition dynamics calculation")

	options = ap.parse_args()

	adata=anndata.read_h5ad(options.input_file)

	# Check for marker genes file and read it into a list
	if bool(options.markers):
		with open(options.markers) as f:
			markergenes = f.read().splitlines()
		markergenes = list(set([sub.replace('-I', '') for sub in markergenes]))

	# Check if user wants to regenerate variable gene selection, or if it needs to be generated from scratch
	if options.hvg == "False":
		if "highly_variable" not in list(adata.var):
			sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	if options.hvg == "True":
			sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	
	# scVelo Core Functions
	scv.pp.filter_and_normalize(adata, min_shared_counts=int(options.minshared), n_top_genes=int(options.topgenes), enforce=True)
	scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
	scv.tl.recover_dynamics(adata, n_jobs=int(options.ncores))
	scv.tl.velocity(adata, mode = 'dynamical')
	scv.tl.velocity_graph(adata)

	# Confirm presence of lower dimensional embeddings and generate if absent
	if options.embedding != "tsne":
		if "X_umap" not in list(adata.obsm):
			scv.tl.umap(adata)
	if options.embedding != "umap":
		if "X_tsne" not in list(adata.obsm):
			scv.tl.tsne(adata)
	scv.tl.louvain(adata)
	scv.pl.velocity_embedding_stream(adata, basis=options.embedding,save="embedding")
	ad.AnnData.write(adata, options.output + "_graph_result.h5ad")

	# Add plotting for Batch Keys if present

	# Check if marker genes are present and plot ones that are
	if bool(options.markers):
		if len(np.setdiff1d(markergenes,adata.var_names)):
			print("Invalid marker genes.")
			print(np.setdiff1d(markergenes,adata.var_names))
			print("were not present in the variable genes list.")
			markergenes = list(set(adata.var_names) & set(markergenes))
		for i in markergenes:
			scv.pl.velocity_embedding_stream(adata, basis=options.embedding, color=[i], save="embedding_"+options.output+"_"+i)

if __name__ == '__main__':
	main()

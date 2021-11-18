import os
import sys
import re
import argparse
import shutil
import warnings
from optparse import OptionParser

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
import matplotlib
import igraph

__author__ = "Anthony S. Castanza"
__email__ = "acastanza@ucsd.edu"
__version__ = "1.0.0"


def main():
    usage = "%prog [options]" + "\n"
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input-file", action="store",
                    dest="input_file", help="Input h5ad (anndata) file from Scanpy.")
    ap.add_argument("-m", "--markers", action="store",
                    dest="markers", nargs='?', help="A list of marker genes")
    ap.add_argument("-v", "--velmode", default="stochastic", action="store", dest="velocity_mode",
                    help="Mode for performing the velocity estimation. 'deterministic', 'stochastic' (default), or 'dynamical'.")
    ap.add_argument("-s", "--shared", default="30", action="store", dest="minshared",
                    help="Filter genes by minimum shared counts.")
    ap.add_argument("-t", "--top", default="2000", action="store", dest="topgenes",
                    help="Top genes for velocity Computation.")
    ap.add_argument("-g", "--hvg", default="seurat_v3", action="store", dest="hvg",
                    help="Recalculate highly_variable_genes.")
    ap.add_argument("-f", "--force", default="False", action="store", dest="enforce",
                    help="Enforce normalizaion using scvelo's internal functions.")
    ap.add_argument("-c", "--pcs", default="30", action="store", dest="pcs",
                    help="Number of principal components used for computing gene moments.")
    ap.add_argument("-n", "--neighbors", default="30", action="store", dest="neighbors",
                    help="Number of nearest neighbors in PCA space used for computing gene moments.")
    ap.add_argument("-d", "--diffkin", default="True", action="store", dest="diff_kinetics",
                    help="Perform differential kinetics analysis using clustering (requires 'dynamical' mode velocity estimation).")
    ap.add_argument("-e", "--embedding", default="umap", action="store", dest="embedding",
                    help="Dataset was processed with 'umap' or 'tsne' embedding.")
    ap.add_argument("-l", "--clustering", default="autodetect_existing", action="store", dest="clustering",
                    help="Some kinetics functions require the dataset to be clustered. Specify 'run_leiden' or 'run_louvain' to create a new clustering, or 'autodetect_existing' to attempt to detect previous clustering with a fallback to 'run_leiden'.")
    ap.add_argument("-r", "--resolution", default="1.0", action="store", dest="resolution",
                    help="Specify a resolution to use for clustering if running the leiden or louvain algorithms.")
    ap.add_argument("-k", "--key_regression", default="NONE", action="store", dest="keys",
                    help="Specify a resolution to use for clustering if running the leiden or louvain algorithms.")
    ap.add_argument("-b", "--batch", default="True", action="store", dest="plot_batches",
                    help="Produce individual velocity plots for each batch in the dataset (if present).")
    ap.add_argument("-p", "--plot", default="png", action="store", dest="plot",
                    help="Save velocity plots as png or svg.")
    ap.add_argument("-o", "--out", default="result", action="store",
                    dest="output", help="Output file basename.")
    ap.add_argument("-j", "--cpu", default="1", action="store", dest="ncores",
                    help="CPU cores to use for transition dynamics calculation.")

    options = ap.parse_args()

    adata = sc.read(options.input_file)

    scv.settings.figdir = "./figures"

    scv.pl.proportions(adata, save=options.output +
                       "_splicing_proportions." + options.plot)

    # Check for marker genes file and read it into a list
    if bool(options.markers):
        with open(options.markers) as f:
            markergenes = f.read().splitlines()
        markergenes = list(set([re.sub('-I$', '', sub)
                                for sub in markergenes]))

    if int(options.topgenes) == 0:
        print("Can't use '0' top genes so resetting to complete gene list. Using " +
              str(len(adata.var)) + " genes.")
        options.topgenes = len(adata.var)

    # Check if user wants to regenerate variable gene selection, or if it needs to be generated from scratch
    if options.hvg == "False":
        if "highly_variable" not in list(adata.var):
            print("Calculation of highly variable genes was not selected but no precomputed set was detected in the dataset so doing it anyway using method " + options.hvg)
            sc.pp.highly_variable_genes(adata, flavor=options.hvg, n_top_genes=int(
                options.topgenes), subset=False, check_values=False)
    else:
        print("Realculation of highly variable genes was requested; calculating using method " + options.hvg)
        sc.pp.highly_variable_genes(adata, flavor=options.hvg, n_top_genes=int(
            options.topgenes), subset=False, check_values=False)

    # scVelo Core Functions
    if options.enforce == "all":
        scv.pp.filter_and_normalize(adata, min_shared_counts=int(
            options.minshared), n_top_genes=int(options.topgenes), subset_highly_variable=False, enforce=True)
    elif options.enforce == "velocity":
        scv.pp.filter_and_normalize(adata, min_shared_counts=int(
            options.minshared), n_top_genes=int(options.topgenes), subset_highly_variable=False, layers_normalize={'spliced', 'unspliced', 'ambiguous'})
    else:
        warnings.warn(
            print("Absolutely no normalization was selected. Using all data layers as-is."))

    if options.keys != "NONE":
        print("Attempting to regress out selected keys: " + options.keys)
        regression_keys = options.keys.split(",")
        sc.pp.regress_out(adata, regression_keys)
        sc.pp.scale(adata)

    if options.hvg == "True" or options.enforce == "all" or options.keys != "NONE":
        print("Additional preprocessing of the non-splicing expression data was requested, re-running PCA and " +
              options.embedding.upper() + " embedding to produce plots that reflect this.")
        sc.pp.pca(adata, n_comps=int(options.pcs), svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=int(options.neighbors))
        if options.embedding == "umap":
            scv.tl.umap(adata)
        if options.embedding == "tsne":
            scv.tl.tsne(adata)

    scv.pp.moments(adata, n_pcs=int(options.pcs),
                   n_neighbors=int(options.neighbors))

    if options.velocity_mode == "dynamical" or options.diff_kinetics == "True":
        scv.tl.recover_dynamics(adata, n_jobs=int(options.ncores))

    scv.tl.velocity(adata, mode=options.velocity_mode)
    scv.tl.velocity_graph(adata, n_jobs=int(options.ncores))
    scv.tl.velocity_pseudotime(adata)

    # Confirm presence of lower dimensional embeddings and generate if absent
    if options.embedding == "umap":
        if "X_umap" not in list(adata.obsm):
            print(
                "'UMAP' Embedding was requested, but we didn't find it in the dataset so creating it now.")
            scv.tl.umap(adata)
    if options.embedding == "tsne":
        if "X_tsne" not in list(adata.obsm):
            print(
                "'tSNE' Embedding was requested, but we didn't find it in the dataset so creating it now.")
            scv.tl.tsne(adata)

    if options.velocity_mode == "dynamical":
        scv.tl.latent_time(adata)

# Detect/create clustering
    if options.clustering == "autodetect_existing":
        if "clusters" in list(adata.obs):
            cluster_type = "clusters"
            cluster_out = "dataset_clusters"
            print(
                "Found 'clusters' key in dataset. We'll use this for plots and any differential kinetics.")
        elif "clusters" not in list(adata.obs):
            if "leiden" in list(adata.obs):
                cluster_type = "leiden"
                cluster_out = "leiden_clusters"
                print(
                    "Found 'leiden' clustering in dataset. We'll use this for plots and any differential kinetics.")
            elif "leiden" not in list(adata.obs):
                if "louvain" in list(adata.obs):
                    cluster_type = "louvain"
                    cluster_out = "louvain_clusters"
                    print(
                        "Found 'louvain' clustering in dataset. We'll use this for plots and any differential kinetics.")
                elif "louvain" not in list(adata.obs):
                    if "walktrap" in list(adata.obs):
                        cluster_type = "walktrap"
                        cluster_out = "walktrap_clusters"
                        print(
                            "Found 'walktrap' clustering in dataset. We'll use this for plots and any differential kinetics.")
                    else:
                        print(
                            "Didn't find any clustering in dataset, clustering data using method: 'leiden'.\nWe'll use this for plots and any differential kinetics.")
                        sc.tl.leiden(adata, resolution=float(
                            options.resolution))
                        cluster_type = "leiden"
                        cluster_out = "leiden_clusters"
    elif options.clustering == "run_louvain":
        print(
            "Clustering data using method: 'louvain'.\nWe'll use this for plots and any differential kinetics.")
        sc.tl.louvain(adata, resolution=float(options.resolution))
        cluster_type = "louvain"
        cluster_out = "louvain_clusters"
    elif options.clustering == "run_leiden":
        print(
            "Clustering data using method: 'leiden'.\nWe'll use this for plots and any differential kinetics.")
        sc.tl.leiden(adata, resolution=float(options.resolution))
        cluster_type = "leiden"
        cluster_out = "leiden_clusters"
    else:
        print(
            "Attempting to use user-specified clustering as-is from key: " + options.clustering)
        cluster_type = options.clustering
        cluster_out = options.clustering + "_clusters"

    sc.tl.rank_genes_groups(adata, cluster_type, method='wilcoxon')

    if options.velocity_mode == "dynamical":
        top_lt_genes = adata.var['fit_likelihood'].sort_values(
            ascending=False).index[:300]
        scv.pl.heatmap(adata, var_names=top_lt_genes, sortby='latent_time', col_color=cluster_type,
                       n_convolve=100, save=options.output + "_top_latent_time_genes_trajectory_heatmap." + options.plot)

    scv.tl.rank_velocity_genes(
        adata, groupby=cluster_type, n_genes=len(adata.var))
    vel_df = scv.DataFrame(adata.uns['rank_velocity_genes']['names']).head(100)
    vel_df.to_csv(options.output + "_top100_velocity_genes_by_" +
                  cluster_out + ".txt", sep="\t")

# # Convert to Gene.By.Sample.Score.Matrix
# unique_values = set()
# for col in scv.DataFrame(adata.uns['rank_velocity_genes']['names']):
#     unique_values.update(scv.DataFrame(
#         adata.uns['rank_velocity_genes']['names'])[col])
# unique_values = list(unique_values)
# unique_values.sort()
#
# gene_by_cluster = pd.DataFrame(columns=scv.DataFrame(
#     adata.uns['rank_velocity_genes']['names']).columns, index=unique_values)
#
# for col in scv.DataFrame(adata.uns['rank_velocity_genes']['names']):
#     gene_by_cluster[col] = list(scv.DataFrame(adata.uns['rank_velocity_genes']['scores'])[
#                    col][np.argsort(scv.DataFrame(adata.uns['rank_velocity_genes']['names'])[col].values)])
#
# gene_by_cluster.to_csv(options.output + "_velocity_gene_scores_by_" +
#                  cluster_out + ".txt", sep="\t")

    if options.velocity_mode == "dynamical":
        scv.tl.rank_dynamical_genes(adata, groupby=cluster_type)
        dyn_df = scv.get_df(adata, 'rank_dynamical_genes/names')
        dyn_df.to_csv(options.output + "_top_dynamical_genes_by_" +
                      cluster_out + ".txt", sep="\t")

# Plotting
    if options.velocity_mode == "dynamical":
        plots = ['latent_time', 'velocity_pseudotime', cluster_type]
    else:
        plots = ['velocity_pseudotime', cluster_type]

    if "batch" in list(adata.obs):
        batches = list(adata.obs['batch'].cat.categories)
        scv.pl.velocity_embedding_stream(
            adata, color=plots + ['batch'], basis=options.embedding, save=options.output + "_velocity_embeddings." + options.plot)
        if options.plot_batches == "True":
            for i in batches:
                try:
                    scv.pl.velocity_embedding_stream(adata[adata.obs['batch'] == i], color=plots, color_map='gnuplot',
                                                     basis=options.embedding, save=options.output + "_batch_" + i + "_velocity_embeddings." + options.plot)
                except ValueError:
                    warnings.warn(print("Unable to plot batch: " + i +
                                        ". Perhaps too many cells were removed by filtering parameters."))
    else:
        scv.pl.velocity_embedding_stream(adata, color=plots, color_map='gnuplot',
                                         basis=options.embedding, save=options.output + "_velocity_embeddings." + options.plot)

    scv.tl.velocity_confidence(adata)
    scv.pl.scatter(adata, c=['velocity_length'], cmap='coolwarm', perc=[
                   5, 95], save=options.output + "_velocity_length_embedding." + options.plot)
    scv.pl.scatter(adata, c=['velocity_confidence'], cmap='coolwarm', perc=[
                   5, 95], save=options.output + "_velocity_confidence_embedding." + options.plot)
    conf_df = adata.obs.groupby(cluster_type)[
        'velocity_length', 'velocity_confidence'].mean().T
    # conf_df.style.background_gradient(cmap='coolwarm', axis=1)
    conf_df.to_csv(options.output + "_velocity_length_and_confidence_by_" +
                   cluster_out + ".txt", sep="\t")

    scv.tl.paga(adata, groups=cluster_type)
    paga_df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
    # paga_df.style.background_gradient(cmap='Blues').format('{:.2g}')
    paga_df.to_csv(options.output + "_paga_transitions_confidence_by_" +
                   cluster_out + ".txt", sep="\t")

    scv.pl.paga(adata, basis=options.embedding, size=50, alpha=.1,
                min_edge_width=2, node_size_scale=1.5, save=options.output + "_paga_velocity_graph_by_" + cluster_out + "." + options.plot)
    scv.pl.scatter(adata, color=['root_cells', 'end_points'],
                   save=options.output + "_velocity_terminal_states." + options.plot)

    ad.AnnData.write(adata, compression="gzip",
                     filename=options.output + "_complete_" + options.velocity_mode + "_velocity_data.h5ad")

    # Check if marker genes are present and plot ones that are
    if options.markers != None:
        scv.settings.figdir = "./figures/markers"
        if len(np.setdiff1d(markergenes, adata.var_names)):
            if "gene_name" in list(adata.var):
                check_names = np.setdiff1d(markergenes, adata.var_names)
                markergenes = list(set(adata.var_names) & set(markergenes))
                for gene in check_names:
                    found_ids = list(
                        adata.var[adata.var["gene_name"] == gene].index)
                    if len(found_ids) > 0:
                        warnings.warn(print("Gene " + gene + " was identified in the dataset under name ",
                                            found_ids, "plots will be created under this ID."))
                        markergenes = markergenes + found_ids
                    else:
                        print("Invalid marker gene: " + gene)
            else:
                print("Invalid marker genes.")
                print(np.setdiff1d(markergenes, adata.var_names))
                print("were not present in the variable genes list.")
                markergenes = list(set(adata.var_names) & set(markergenes))
        for i in markergenes:
            scv.pl.velocity_embedding_stream(adata, basis=options.embedding, color=[
                i], save=options.output + "_embedding_" + i + "." + options.plot)
        scv.pl.velocity(adata, markergenes, ncols=1,
                        save=options.output + "_combined_per-marker_velocity." + options.plot)

   # Stuff for Differential Kinetics
    if options.diff_kinetics == "True":
        scv.settings.figdir = "./diff_kinetics_figures"
        velocity_genes_list = list(
            adata.var['velocity_genes'][adata.var['velocity_genes'] == True].index)
        scv.tl.differential_kinetic_test(adata, groupby=cluster_type)
        kdf = scv.get_df(adata[:, velocity_genes_list], [
                         'fit_diff_kinetics', 'fit_pval_kinetics'], precision=2)
        kdf.to_csv(options.output + "_differential_kinetics_test_for_velocity_genes_by_" +
                   cluster_out + ".txt", sep="\t")
        kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
        top_genes = adata.var['fit_likelihood'].sort_values(
            ascending=False).index[:100]
        scv.pl.scatter(adata, basis=top_genes[:20], ncols=5, add_outline='fit_diff_kinetics', **kwargs,
                       save=options.output + "_top20_fit_likelihood_genes_after_differential_kinetics." + options.plot)
        diff_clusters = list(
            adata[:, velocity_genes_list].var['fit_diff_kinetics'])
        scv.pl.scatter(adata, size=60, title='diff kinetics',
                       add_outline=diff_clusters, outline_width=(.8, .2), color=cluster_type, save=options.output + "_outlined_clusters_affected_by_differential_kinetics." + options.plot)
        scv.tl.velocity(adata, mode=options.velocity_mode, diff_kinetics=True)
        scv.tl.velocity_graph(adata, n_jobs=int(options.ncores))
        scv.tl.velocity_pseudotime(adata)
        if options.velocity_mode == "dynamical":
            scv.tl.latent_time(adata)
            top_lt_genes = adata.var['fit_likelihood'].sort_values(
                ascending=False).index[:300]
            scv.pl.heatmap(adata, var_names=top_lt_genes, sortby='latent_time', col_color=cluster_type, n_convolve=100,
                           save=options.output + "_top_latent_time_genes_trajectory_heatmap_after_differential_kinetics." + options.plot)

        scv.tl.rank_velocity_genes(
            adata, groupby=cluster_type, n_genes=len(adata.var))
        vel_dk_df = scv.DataFrame(
            adata.uns['rank_velocity_genes']['names']).head(100)
        vel_dk_df.to_csv(options.output + "_top100_velocity_genes_by_" +
                         cluster_out + "_after_differential_kinetics.txt", sep="\t")

        if options.velocity_mode == "dynamical":
            scv.tl.rank_dynamical_genes(adata, groupby=cluster_type)
            dyn_dk_df = scv.get_df(adata, 'rank_dynamical_genes/names')
            dyn_dk_df.to_csv(options.output + "_top_dynamical_genes_by_" +
                             cluster_out + "_after_differential_kinetics.txt", sep="\t")

        if "batch" in list(adata.obs):
            batches = list(adata.obs['batch'].cat.categories)
            scv.pl.velocity_embedding_stream(
                adata, color=plots + ['batch'], basis=options.embedding, save=options.output + "_velocity_embeddings_after_differential_kinetics." + options.plot)
            if options.plot_batches == "True":
                for i in batches:
                    try:
                        scv.pl.velocity_embedding_stream(adata[adata.obs['batch'] == i], color=plots, color_map='gnuplot',
                                                         basis=options.embedding, save=options.output + i + "_velocity_embeddings_after_differential_kinetics." + options.plot)
                    except ValueError:
                        warnings.warn(print("Unable to plot batch: " + i +
                                            ". Perhaps too many cells were removed by filtering parameters."))
        else:
            scv.pl.velocity_embedding_stream(adata, color=plots, color_map='gnuplot',
                                             basis=options.embedding, save=options.output + "_velocity_embeddings_after_differential_kinetics." + options.plot)

        scv.tl.velocity_confidence(adata)
        scv.pl.scatter(adata, c=['velocity_length'], cmap='coolwarm', perc=[
                       5, 95], save=options.output + "_velocity_length_embedding_after_differential_kinetics." + options.plot)
        scv.pl.scatter(adata, c=['velocity_confidence'], cmap='coolwarm', perc=[
                       5, 95], save=options.output + "_velocity_confidence_embedding_after_differential_kinetics." + options.plot)
        conf_dk_df = adata.obs.groupby(cluster_type)[
            'velocity_length', 'velocity_confidence'].mean().T
        # conf_dk_df.style.background_gradient(cmap='coolwarm', axis=1)
        conf_dk_df.to_csv(options.output + "_velocity_length_and_confidence_by_" +
                          cluster_out + "_after_differential_kinetics.txt", sep="\t")

        scv.tl.paga(adata, groups=cluster_type)
        paga_dk_df = scv.get_df(
            adata, 'paga/transitions_confidence', precision=2).T
        # paga_dk_df.style.background_gradient(cmap='Blues').format('{:.2g}')
        paga_dk_df.to_csv(options.output + "_paga_transitions_confidence_by_" +
                          cluster_out + "_after_differential_kinetics.txt", sep="\t")

        scv.pl.paga(adata, basis=options.embedding, size=50, alpha=.1,
                    min_edge_width=2, node_size_scale=1.5, save=options.output + "_paga_velocity_graph_by_" + cluster_out + "_after_differential_kinetics." + options.plot)
        scv.pl.scatter(adata, color=['root_cells', 'end_points'], save=options.output +
                       "_velocity_terminal_states_after_differential_kinetics." + options.plot)

        ad.AnnData.write(adata, compression="gzip",
                         filename=options.output + "_complete_" + options.velocity_mode + "_velocity_post-differential_kinetics_data.h5ad")

        if options.markers != None:
            scv.settings.figdir = "./diff_kinetics_figures/markers"
            for i in markergenes:
                scv.pl.velocity_embedding_stream(adata, basis=options.embedding, color=[
                    i], save="Marker_" + i + "_" + options.output + "_embedding_after_differential_kinetics." + options.plot)
            scv.pl.velocity(adata, markergenes, ncols=1,
                            save=options.output + "_combined_per-marker_velocity_after_differential_kinetics." + options.plot)


if __name__ == '__main__':
    main()

# Functions to Enable Downstream Analysis of Gene Sets for Single Cell Datasets

# Commands to easily work with these functions
# import sys, importlib
# sys.path.insert(1, '/Users/acastanza/github/scVelo.ComputeVelocity/module')
# import anndata as ad
# import GeneSetAnalysisFunctions
# importlib.reload(GeneSetAnalysisFunctions)
# adata = ad.read_h5ad(options.input_file)
# import matplotlib
# matplotlib.use('Agg')
import sys, math, re, numpy, pandas, scvelo
from scipy.sparse import isspmatrix
import GeneSetAnalysisFunctions


# Convert to Gene.By.Sample.Score.Matrix
def get_gene_values(adata, key='X', genes_min_nonzero_cells=0, outname="Dataset", velocity_weight=False, write_gct=True):
    if key.upper() == 'X':
        use_cell_level_filter = True
        if isspmatrix(adata.X):
            gene_by_cell = pandas.DataFrame(
                adata.X.todense()).transpose()
        else:
            gene_by_cell = pandas.DataFrame(
                adata.X).transpose()
        gene_by_cell.index = adata.var.index
        gene_by_cell.columns = adata.obs.index
        out_matrix = gene_by_cell
        filename = outname + "_" + "cell_level_genes_" + key
    elif key in adata.layers:
        use_cell_level_filter = True
        if isspmatrix(adata.layers[key]):
            gene_by_cell = pandas.DataFrame(
                adata.layers[key].todense()).transpose()
        else:
            gene_by_cell = pandas.DataFrame(adata.layers[key].transpose())
        gene_by_cell.index = adata.var.index
        gene_by_cell.columns = adata.obs.index
        out_matrix = gene_by_cell
        filename = outname + "_" + "cell_level_genes_" + key
    elif key == 'rank_velocity_genes' or key == 'rank_genes_groups':
        # key='rank_velocity_genes' and key='rank_genes_groups' both work
        use_cell_level_filter = False
        cluster_key = detect_clusters(adata)
        unique_values = set()
        for col in scvelo.DataFrame(adata.uns[key]['names']):
            unique_values.update(scvelo.DataFrame(
                adata.uns[key]['names'])[col])
        unique_values = list(unique_values)
        unique_values.sort()
        gene_by_cluster = pandas.DataFrame(columns=scvelo.DataFrame(
            adata.uns[key]['names']).columns, index=unique_values)
        for col in scvelo.DataFrame(adata.uns[key]['names']):
            gene_by_cluster[col] = pandas.DataFrame({'scores':scvelo.DataFrame(adata.uns[key]['scores'])[col].values},index=scvelo.DataFrame(adata.uns[key]['names'])[col].values).reindex(gene_by_cluster.index).fillna(0)
        rank_genes_groups_by_cluster = gene_by_cluster.copy()
        if velocity_weight == True and key != 'rank_velocity_genes':
            if "rank_velocity_genes" in adata.uns:
                print("Applying velocity weights to ranked genes")
                unique_values = set()
                for col in scvelo.DataFrame(adata.uns['rank_velocity_genes']['names']):
                    unique_values.update(scvelo.DataFrame(
                    adata.uns['rank_velocity_genes']['names'])[col])
                unique_values = list(unique_values)
                unique_values.sort()
                velocity_gene_by_cluster = pandas.DataFrame(columns=scvelo.DataFrame(
                    adata.uns['rank_velocity_genes']['names']).columns, index=unique_values)
                for col in scvelo.DataFrame(adata.uns['rank_velocity_genes']['names']):
                    velocity_gene_by_cluster[col] = pandas.DataFrame({'scores':scvelo.DataFrame(adata.uns['rank_velocity_genes']['scores'])[col].values},index=scvelo.DataFrame(adata.uns['rank_velocity_genes']['names'])[col].values).reindex(velocity_gene_by_cluster.index).fillna(0)
                rank_velocity_genes_by_cluster = velocity_gene_by_cluster.copy()
                velocity_weight = (1 + numpy.log(1 + numpy.absolute(rank_velocity_genes_by_cluster.reindex(rank_genes_groups_by_cluster.index).fillna(0))))
                out_matrix = numpy.sign(rank_genes_groups_by_cluster) * (numpy.absolute(rank_genes_groups_by_cluster) ** velocity_weight)
                filename = outname + "_" + "velocity_weighted_ranked_genes"
        else:
            out_matrix = rank_genes_groups_by_cluster
            filename = outname + "_" + key + "_by_" + cluster_key + "_clusters"
    out_matrix.index.name = "NAME"
    out_matrix.index = out_matrix.index.str.replace(
        '\\..*', '', regex=True)
    if use_cell_level_filter == True:
        if int(genes_min_nonzero_cells) > 0:
            out_matrix = out_matrix[out_matrix.mask(out_matrix!=0).count(axis=1) > int(genes_min_nonzero_cells)]
    out_matrix.insert(loc=0, column='Description', value="NA")
    if write_gct == True:
        text_file = open(filename + ".gct", "w")
        text_file.write('#1.2\n')
        text_file.write(str(len(out_matrix)) + "\t" +
                        str(len(out_matrix.columns) - 1) + "\n")
        text_file.close()
        out_matrix.to_csv(filename + ".gct", sep="\t", mode='a')
    return {'data': out_matrix.drop(labels="Description", axis=1), 'row_descriptions': out_matrix["Description"].values, 'outname': filename}
    # sumtest=Dataset_rank_velocity_genes.reindex(Dataset_rank_genes_groups.index).fillna(0) + Dataset_rank_genes_groups


def detect_clusters(adata, silent=True):
    if "clusters" in list(adata.obs):
        cluster_type = "clusters"
        if silent == False:
            print(
                "Found 'clusters' key in dataset")
    elif "clusters" not in list(adata.obs):
        if "leiden" in list(adata.obs):
            cluster_type = "leiden"
            if silent == False:
                print(
                    "Found 'leiden' clustering in dataset.")
        elif "leiden" not in list(adata.obs):
            if "louvain" in list(adata.obs):
                cluster_type = "louvain"
                if silent == False:
                    print(
                        "Found 'louvain' clustering in dataset.")
            elif "louvain" not in list(adata.obs):
                if "walktrap" in list(adata.obs):
                    cluster_type = "walktrap"
                    if silent == False:
                        print(
                            "Found 'walktrap' clustering in dataset.")
                else:
                    print("No clustering found in the dataset.")
                    sys.exit(1)
    return cluster_type


def make_pseudobulk(adata, key='X', method="sum", genes_min_nonzero_cells=0, clustering="detect", outname="Dataset", velocity_weight=False, write_gct=True):
    gene_values = get_gene_values(adata, key='X', write_gct=False)
    gene_values = gene_values['data']
    if int(genes_min_nonzero_cells) > 0:
        gene_values = gene_values[gene_values.mask(gene_values!=0).count(axis=1) > int(genes_min_nonzero_cells)]
    if clustering == "detect":
        clustering = detect_clusters(adata, silent=True)
    cluster_assignments = adata.obs[clustering]
    gene_values = gene_values.transpose()
    gene_values.insert(loc=0, column='Clusters', value=cluster_assignments)
    if method == "sum":
        pseudobulk_df=gene_values.groupby(["Clusters"]).sum()
    if method == "mean":
        pseudobulk_df=gene_values.groupby(["Clusters"]).mean()
    if method == "median":
        pseudobulk_df=gene_values.groupby(["Clusters"]).median()
    if method == "max":
        pseudobulk_df=gene_values.groupby(["Clusters"]).max()
    pseudobulk_df = pseudobulk_df.transpose()
    pseudobulk_df.columns=pseudobulk_df.columns.to_list()
    if velocity_weight == True:
        if "rank_velocity_genes" in adata.uns:
            print("Applying velocity weights to pseudobulk counts")
            unique_values = set()
            for col in scvelo.DataFrame(adata.uns['rank_velocity_genes']['names']):
                unique_values.update(scvelo.DataFrame(
                    adata.uns['rank_velocity_genes']['names'])[col])
            unique_values = list(unique_values)
            unique_values.sort()
            velocity_gene_by_cluster = pandas.DataFrame(columns=scvelo.DataFrame(
                adata.uns['rank_velocity_genes']['names']).columns, index=unique_values)
            for col in scvelo.DataFrame(adata.uns['rank_velocity_genes']['names']):
                velocity_gene_by_cluster[col] = pandas.DataFrame({'scores':scvelo.DataFrame(adata.uns['rank_velocity_genes']['scores'])[col].values},index=scvelo.DataFrame(adata.uns['rank_velocity_genes']['names'])[col].values).reindex(velocity_gene_by_cluster.index).fillna(0)
            rank_velocity_genes_by_cluster = velocity_gene_by_cluster.copy()
            velocity_weight = (1 + numpy.log(1 + numpy.absolute(rank_velocity_genes_by_cluster.reindex(pseudobulk_df.index).fillna(0))))
            out_matrix = numpy.sign(pseudobulk_df) * (numpy.absolute(pseudobulk_df) ** velocity_weight)
            filename = outname + "_" + "_cluster_level_velocity_weighted_pseudobulk_counts"
    else:
        out_matrix = pseudobulk_df
        filename = outname + "_cluster_level_pseudobulk_counts"
    out_matrix.index.name="NAME"
    out_matrix.insert(loc=0, column='Description', value="NA")
    if write_gct == True:
        text_file = open(filename + ".gct", "w")
        text_file.write('#1.2\n')
        text_file.write(str(len(out_matrix)) + "\t" +
                        str(len(out_matrix.columns) - 1) + "\n")
        text_file.close()
        out_matrix.to_csv(filename + ".gct", sep="\t", mode='a')
    return {'data': out_matrix.drop(labels="Description", axis=1), 'row_descriptions': out_matrix["Description"].values, 'outname': filename}


def load_ssgsea_result(ssgsea_result):
    ssgsea_df = pandas.read_csv(ssgsea_result, sep='\t', header=2, index_col=[
                            0, 1], skip_blank_lines=True)
    ssgsea_df.index = ssgsea_df.index.droplevel(1)  # Drop gene descriptions
    return ssgsea_df


# Add Clusterwise ssGSEA scores to adata.obs as a cell level score for plotting
def expand_ssgsea_cluster_scores(adata, cluster_key, ssgsea_result):
    ssgsea_df = load_ssgsea_result(ssgsea_result)
    ssgsea_cell_df = ssgsea_df.transpose()
    ssgsea_cell_df = ssgsea_cell_df.reindex(list(adata.obs[cluster_key]))
    ssgsea_cell_df.index = range(len(ssgsea_cell_df.index))
    ssgsea_cell_df = ssgsea_cell_df.set_index(adata.obs.index)
    adata.uns['gene_sets']=list(ssgsea_cell_df.columns)
    adata.obs[ssgsea_cell_df.columns] = ssgsea_cell_df[ssgsea_cell_df.columns]
    return ssgsea_cell_df


def import_ssgsea_cell_scores(adata, ssgsea_result):
    ssgsea_df = load_ssgsea_result(ssgsea_result)
    ssgsea_cell_df = ssgsea_df.transpose()
    ssgsea_cell_df = ssgsea_cell_df.reindex(list(adata.obs.index))
    adata.uns['gene_sets']=list(ssgsea_cell_df.columns)
    adata.obs[ssgsea_cell_df.columns] = ssgsea_cell_df[ssgsea_cell_df.columns]
    return ssgsea_cell_df


def ssgsea_plot_all(adata, ssgsea_result, basis, outname, format):  # Plotting
    cluster_key = detect_clusters(adata)
    ssgsea_cell_df = expand_ssgsea_cluster_scores(
        adata, cluster_key, ssgsea_result)
    ssgsea_sets = list(ssgsea_cell_df.columns)
    for set in ssgsea_sets:
        scvelo.pl.velocity_embedding_stream(adata, basis=basis, color=[
                                         set, cluster_key], color_map='seismic', save=set + "_" + outname + "_embedding." + format)


def ssgsea_plot_hits(adata, filtered_set_hits, ssgsea_result, basis, outname="dataset", format="png"):  # Plotting
    cluster_key = detect_clusters(adata)
    ssgsea_cell_df = expand_ssgsea_cluster_scores(
        adata, cluster_key, ssgsea_result)
    for i in range(len(filtered_set_hits)):
        set = str(filtered_set_hits.index[i])
        scvelo.pl.velocity_embedding_stream(adata, basis=basis, color=[
            set, cluster_key], color_map='seismic', add_outline=[filtered_set_hits.iloc[i][0], filtered_set_hits.iloc[i][2]], save=outname + "_" + set + "_Cluster_" + str(filtered_set_hits.iloc[i][0]) + "_to_" + str(filtered_set_hits.iloc[i][2]) + "_" + basis + "_embedding." + format)


# Calculate the Gene Set ES Delta pairwise for every set
def create_transition_matrix(ssgsea_result, set):
    ssgsea_raw_df = load_ssgsea_result(ssgsea_result)
    ssgsea_sets = list(ssgsea_raw_df.index)
    set_transition = pandas.DataFrame(
        columns=ssgsea_raw_df.columns, index=ssgsea_raw_df.columns)
    if len(numpy.where(ssgsea_raw_df.index == set)) == 1:
        test_set = ssgsea_raw_df.iloc[[
            int(numpy.where(ssgsea_raw_df.index == set)[0])]]
    else:
        print("Found Duplicate Gene Set Name: " + set)
        sys.exit(1)
    for first_cluster in test_set.columns:
        for second_cluster in test_set.columns:
            set_transition.at[first_cluster, second_cluster] = float(
                test_set[second_cluster]) - float(test_set[first_cluster])
    return set_transition


# Take the pairwise deltas from create_transition_matrix and screen out invalid transitions using the PAGA matrix as a mask then apply thresholding criteria
def find_candidate_transitions(adata, ssgsea_result, set, conf_threshold=0.3, adj_threshold=0.1, stdev_filter=[2], silent=False):
    # connectivities confidence
    paga_conf_df = scvelo.get_df(
        adata, 'paga/transitions_confidence', precision=4).T
    # connectivities adjacency
    paga_adj_df = scvelo.get_df(adata, 'paga/connectivities', precision=4).T
    # paga_tree_df = scvelo.get_df(adata, 'paga/connectivities_tree', precision=2).T # connectivities subtree
    cluster_key = detect_clusters(adata)
    if set not in adata.obs.columns:
        ssgsea_cell_df = expand_ssgsea_cluster_scores(
            adata, cluster_key, ssgsea_result)
    set_transition = create_transition_matrix(
        ssgsea_result, set)
    set_transition_full = set_transition.copy()
    set_transition_full[:] = numpy.where(numpy.arange(set_transition_full.shape[0])[
                                      :, None] >= numpy.arange(set_transition_full.shape[1]), numpy.nan, set_transition_full)
    # Maybe apply the paga_adj_df here to make sure the standard deviations are based on just the adjacent clusters
    set_transition_full_list = set_transition_full.values.tolist()
    flat_set_transition_full_list = [
        item for sublist in set_transition_full_list for item in sublist if math.isnan(item) == False]
    flat_set_transition_full_list = list(
        map(abs, flat_set_transition_full_list))
    set_transition_pass = set_transition[paga_conf_df[paga_adj_df > float(
        adj_threshold)] > float(conf_threshold)]
    set_transition_pass_abs = set_transition_pass.abs().copy()
    mean = numpy.mean(flat_set_transition_full_list)
    standard_deviation = numpy.std(flat_set_transition_full_list)
    distance_from_mean = abs(set_transition_pass_abs - mean)
    if len(stdev_filter) == 2:
        filtered = numpy.logical_and(distance_from_mean < (
            max(stdev_filter) * standard_deviation), distance_from_mean > (min(stdev_filter) * standard_deviation))
    if len(stdev_filter) == 1:
        filtered = distance_from_mean > (2 * standard_deviation)
    filtered_locs = list(numpy.where(filtered))
    transition_locs = list(filtered[filtered == True].stack().index)
    ssgsea_raw_df = load_ssgsea_result(ssgsea_result)
    test_set = ssgsea_raw_df.iloc[[
        int(numpy.where(ssgsea_raw_df.index == set)[0])]]
    set_hits = []
    for i in range(len(transition_locs)):
        if silent == False:
            print("Gene set " + set + " was scored as a candidate for transition from Cluster " + str(transition_locs[i][0]) + " (Enrichment Score: " + str(test_set.loc[set, transition_locs[i][0]].round(2)) + ") to Cluster " + str(
                transition_locs[i][1]) + " (Enrichment Score: " + str(test_set.loc[set, transition_locs[i][1]].round(2)) + ") at PAGA transition confidence >" + str(conf_threshold) + " and adjacency >" + str(adj_threshold))
        set_hits.append([set, str(transition_locs[i][0]), test_set.loc[set, transition_locs[i][0]], str(transition_locs[i][1]),
                         test_set.loc[set, transition_locs[i][1]], test_set.loc[set, transition_locs[i][1]] - test_set.loc[set, transition_locs[i][0]]])
    return set_hits


# Using the results of find_candidate_transitions keep candidates that have good directionality
def find_good_transitions(adata, ssgsea_result, conf_threshold=0.3, adj_threshold=0.1, stdev_filter=[2], silent=False):
    ssgsea_raw_df = load_ssgsea_result(ssgsea_result)
    all_sets = ssgsea_raw_df.index.to_list()
    all_set_results = []
    for set in all_sets:
        set_hits = find_candidate_transitions(
            adata=adata, ssgsea_result=ssgsea_result, set=set, conf_threshold=conf_threshold, adj_threshold=adj_threshold, stdev_filter=stdev_filter, silent=silent)
        all_set_results.append(set_hits)
    all_set_results_flat = [
        item for sublist in all_set_results for item in sublist]
    all_set_results_df = pandas.DataFrame(all_set_results_flat)
    # Sets have a Positive Change
    all_positive_changes = all_set_results_df[all_set_results_df[5] > 0]
    all_positive_changes = all_positive_changes[all_positive_changes[4].astype(
        float) > 0]  # Ending Cluster Ends Positive
    all_positive_changes.columns = ["Gene_Set", "Start_Cluster",
                                    "Start_Cluster_ES", "End_Cluster", "End_Cluster_ES", "Cluster_ES_Delta"]
    all_positive_changes.sort_values(
        by="Cluster_ES_Delta", ascending=False, inplace=True)
    # Sets have a Negative Change
    all_negative_changes = all_set_results_df[all_set_results_df[5] < 0]
    all_negative_changes = all_negative_changes[all_negative_changes[2].astype(
        float) > 0]  # Starting Cluster Starts Positive
    all_negative_changes.columns = ["Gene_Set", "Start_Cluster",
                                    "Start_Cluster_ES", "End_Cluster", "End_Cluster_ES", "Cluster_ES_Delta"]
    all_negative_changes.sort_values(
        by="Cluster_ES_Delta", ascending=False, inplace=True)
    filtered_set_hits = all_positive_changes.append(all_negative_changes)
    filtered_set_hits.set_index("Gene_Set", inplace=True)
    return filtered_set_hits


def str_to_bool(s):
    if s == 'True':
         return True
    elif s == True:
        return True
    elif s == 'False':
             return False
    elif s == False:
        return False
    else:
         raise ValueError # evil ValueError that doesn't tell you what the wrong value was


# Get the two clusters latent time and cor() latent time with the clusters ES's. for PAGA transitions, (then compute threholds?)
# import nympy as np
# import scipy
# time_and_cluster_per_cell = adata.obs[["leiden","velocity_pseudotime"]]
# clusters = list(set(time_and_cluster_per_cell["leiden"]))
# time_per_cluster = []
# for cluster in clusters:
#     time_per_cluster.append(time_and_cluster_per_cell["velocity_pseudotime"][time_and_cluster_per_cell["leiden"]==cluster].median())
# cluster_times = dict(zip(clusters,time_per_cluster))
# numpy.asarray([float(cluster_times.get("0")),float(cluster_times.get("1"))])
# for i in range(len(filtered_set_hits)):
#     set = str(filtered_set_hits.index[i])
#  scipy.stats.pearsonr(numpy.asarray([float(cluster_times.get("0")),float(cluster_times.get("1"))]),numpy.array([float(filtered_set_hits.iat[0,1]),float(filtered_set_hits.iat[0,3])]))


# GeneSetAnalysisFunctions.ssgsea_plot_hits(adata,GeneSetAnalysisFunctions.find_good_transitions(adata,"/Users/acastanza/Downloads/E14_5_Pancreas_dim_reduce_clustered_complete_stochastic_velocity_data_velocity_weighted_ranked_genes.PROJ.gct", conf_threshold=0.2, adj_threshold=0.1, stdev_filter=[1.5]),"/Users/acastanza/Downloads/E14_5_Pancreas_dim_reduce_clustered_complete_stochastic_velocity_data_velocity_weighted_ranked_genes.PROJ.gct", basis="umap", outname="E14_5_Pancreas_velocity_weighted_enrichment")

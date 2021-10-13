# Functions to Enable Downstream Analysis of Gene Sets for Single Cell Datasets

# Commands to easily work with these functions
# import sys, importlib
# sys.path.insert(1, '/Users/acastanza/github/scVelo.ComputeVelocity/module')
# import anndata as ad
# import GeneSetAnalysisFunctions
# importlib.reload(GeneSetAnalysisFunctions)
# adata = ad.read_h5ad(options.input_file)

# Convert to Gene.By.Sample.Score.Matrix
def velocity_score_to_gct(adata, outkey='ranked_velocity_genes', outname="Dataset"):
    import re
    import numpy as np
    import pandas as pd
    import GeneSetAnalysisFunctions
    import scvelo as scv
    unique_values = set()
    for col in scv.DataFrame(adata.uns[outkey]['names']):
        unique_values.update(scv.DataFrame(
            adata.uns[outkey]['names'])[col])
    unique_values = list(unique_values)
    unique_values.sort()
    gene_by_cluster = pd.DataFrame(columns=scv.DataFrame(
        adata.uns[outkey]['names']).columns, index=unique_values)
    for col in scv.DataFrame(adata.uns[outkey]['names']):
        gene_by_cluster[col] = list(scv.DataFrame(adata.uns[outkey]['scores'])[
            col][np.argsort(scv.DataFrame(adata.uns[outkey]['names'])[col].values)])
    gene_by_cluster.index.name = "NAME"
    gene_by_cluster.index = gene_by_cluster.index.str.replace(
        '\\..*', '', regex=True)
    cluster_key = GeneSetAnalysisFunctions.detect_clusters(adata)
    gene_by_cluster.insert(loc=0, column='Description', value="NA")
    text_file = open(outname + "_" + outkey + "_by_" +
                     cluster_key + "_clusters.gct", "w")
    text_file.write('#1.2\n')
    text_file.write(str(len(gene_by_cluster)) + "\t" +
                    str(len(gene_by_cluster.columns) - 1) + "\n")
    text_file.close()
    gene_by_cluster.to_csv(outname + "_" + outkey + "_by_" +
                           cluster_key + "_clusters.gct", sep="\t", mode='a')


def load_ssgsea_result(ssgsea_result):
    import sys
    import pandas as pd
    ssgsea_df = pd.read_csv(ssgsea_result, sep='\t', header=2, index_col=[
                            0, 1], skip_blank_lines=True)
    ssgsea_df.index = ssgsea_df.index.droplevel(1)  # Drop gene descriptions
    return ssgsea_df


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


# Add Clusterwise ssGSEA scores to adata.obs as a cell level score for plotting
def adata_import_ssgsea_scores(adata, cluster_key, ssgsea_result):
    import GeneSetAnalysisFunctions
    ssgsea_df = GeneSetAnalysisFunctions.load_ssgsea_result(ssgsea_result)
    ssgsea_cell_df = ssgsea_df.transpose()
    ssgsea_cell_df = ssgsea_cell_df.reindex(list(adata.obs[cluster_key]))
    ssgsea_cell_df.index = range(len(ssgsea_cell_df.index))
    ssgsea_cell_df = ssgsea_cell_df.set_index(adata.obs.index)
    adata.obs[ssgsea_cell_df.columns] = ssgsea_cell_df[ssgsea_cell_df.columns]
    return ssgsea_cell_df


def ssgsea_plot_all(adata, ssgsea_result, basis, outname, format):  # Plotting
    import GeneSetAnalysisFunctions
    import scvelo as scv
    cluster_key = GeneSetAnalysisFunctions.detect_clusters(adata)
    ssgsea_cell_df = GeneSetAnalysisFunctions.adata_import_ssgsea_scores(
        adata, cluster_key, ssgsea_result)
    ssgsea_sets = list(ssgsea_cell_df.columns)
    for set in ssgsea_sets:
        scv.pl.velocity_embedding_stream(adata, basis=basis, color=[
                                         set, cluster_key], color_map='seismic', save=set + "_" + outname + "_embedding." + format)


def ssgsea_plot_hits(adata, filtered_set_hits, ssgsea_result, basis, outname, format):  # Plotting
    import GeneSetAnalysisFunctions
    import scvelo as scv
    cluster_key = GeneSetAnalysisFunctions.detect_clusters(adata)
    ssgsea_cell_df = GeneSetAnalysisFunctions.adata_import_ssgsea_scores(
        adata, cluster_key, ssgsea_result)
    for i in range(len(filtered_set_hits)):
        set = str(filtered_set_hits.index[i])
        scv.pl.velocity_embedding_stream(adata, basis=basis, color=[
            set, cluster_key], color_map='seismic', add_outline=[filtered_set_hits.iloc[i][0], filtered_set_hits.iloc[i][2]], save=set + "_Cluster_" + str(filtered_set_hits.iloc[i][0]) + "_to_" + str(filtered_set_hits.iloc[i][2]) + "_" + outname + "_embedding." + format)


# Calculate the Gene Set ES Delta pairwise for every set
def create_transition_matrix(ssgsea_result, set):
    import GeneSetAnalysisFunctions
    import sys
    import pandas as pd
    import numpy as np
    ssgsea_raw_df = GeneSetAnalysisFunctions.load_ssgsea_result(ssgsea_result)
    ssgsea_sets = list(ssgsea_raw_df.index)
    set_transition = pd.DataFrame(
        columns=ssgsea_raw_df.columns, index=ssgsea_raw_df.columns)
    if len(np.where(ssgsea_raw_df.index == set)) == 1:
        test_set = ssgsea_raw_df.iloc[[
            int(np.where(ssgsea_raw_df.index == set)[0])]]
    else:
        print("Found Duplicate Gene Set Name: " + set)
        sys.exit(1)
    for first_cluster in test_set.columns:
        for second_cluster in test_set.columns:
            set_transition.at[first_cluster, second_cluster] = float(
                test_set[second_cluster]) - float(test_set[first_cluster])
    return set_transition


# Take the pairwise deltas from create_transition_matrix and screen out invalid transitions using the PAGA matrix as a mask then apply thresholding criteria
def find_candidate_transitions(adata, ssgsea_result, set, conf_threshold=0.5, adj_threshold=0.5, silent=False):
    import GeneSetAnalysisFunctions
    import numpy as np
    import pandas as pd
    import scvelo as scv
    import math
    # connectivities confidence
    paga_conf_df = scv.get_df(
        adata, 'paga/transitions_confidence', precision=2).T
    # connectivities adjacency
    paga_adj_df = scv.get_df(adata, 'paga/connectivities', precision=2).T
    # paga_tree_df = scv.get_df(adata, 'paga/connectivities_tree', precision=2).T # connectivities subtree
    cluster_key = GeneSetAnalysisFunctions.detect_clusters(adata)
    ssgsea_cell_df = GeneSetAnalysisFunctions.adata_import_ssgsea_scores(
        adata, cluster_key, ssgsea_result)
    set_transition = GeneSetAnalysisFunctions.create_transition_matrix(
        ssgsea_result, set)
    set_transition_pass = set_transition[paga_conf_df[paga_adj_df > float(
        adj_threshold)] > float(conf_threshold)]
    set_transition_pass_list = set_transition_pass.values.tolist()
    flat_set_transition_pass_list = [
        item for sublist in set_transition_pass_list for item in sublist if math.isnan(item) == False]
    mean = np.mean(flat_set_transition_pass_list)
    standard_deviation = np.std(flat_set_transition_pass_list)
    distance_from_mean = abs(flat_set_transition_pass_list - mean)
    max_deviations = 2
    filtered = np.logical_and(distance_from_mean < (
        max_deviations * standard_deviation), distance_from_mean > (1 * standard_deviation))
    filtered_locs = list(np.where(filtered)[0])
    transition_values = np.array(
        flat_set_transition_pass_list)[filtered_locs]
    transition_values = list(transition_values)
    ssgsea_raw_df = GeneSetAnalysisFunctions.load_ssgsea_result(ssgsea_result)
    test_set = ssgsea_raw_df.iloc[[
        int(np.where(ssgsea_raw_df.index == set)[0])]]
    set_hits = []
    for i in range(len(transition_values)):
        if silent == False:
            print("Gene set " + set + " was scored as a candidate for transition from Cluster " + str(np.where(set_transition_pass == transition_values[i])[0]) + " (Enrichment Score: " + str(float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[0]))])) + ") to Cluster " + str(np.where(
                set_transition_pass == transition_values[i])[1]) + " (Enrichment Score: " + str(float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[1]))])) + ") at PAGA transition confidence >" + str(conf_threshold) + " and adjacency >" + str(adj_threshold))
        set_hits.append([set, str(np.where(set_transition_pass == transition_values[i])[0]).strip("[]"), str(float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[0]))])), str(np.where(set_transition_pass == transition_values[i])[1]).strip("[]"), str(
            float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[1]))])), float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[1]))]) - float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[0]))])])
    return(set_hits)


# Using the results of find_candidate_transitions keep candidates that have good directionality
def find_good_transitions(adata, ssgsea_result, conf_threshold=0.5, adj_threshold=0.5, silent=False):
    import GeneSetAnalysisFunctions
    import pandas as pd
    ssgsea_raw_df = GeneSetAnalysisFunctions.load_ssgsea_result(ssgsea_result)
    all_sets = ssgsea_raw_df.index.to_list()
    all_set_results = []
    for set in all_sets:
        set_hits = GeneSetAnalysisFunctions.find_candidate_transitions(
            adata=adata, ssgsea_result=ssgsea_result, set=set, conf_threshold=0.5, adj_threshold=0.5, silent=silent)
        all_set_results.append(set_hits)
    all_set_results_flat = [
        item for sublist in all_set_results for item in sublist]
    all_set_results_df = pd.DataFrame(all_set_results_flat)
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
    return(filtered_set_hits)

# Get the two clusters latent time and cor() latent time with the clusters ES's. for PAGA transitions, (then compute threholds?)
# time_and_cluster_per_cell = adata.obs[["leiden","velocity_pseudotime"]]
# clusters = list(set(time_and_cluster_per_cell["leiden"]))
# time_per_cluster = []
# for cluster in clusters:
#     time_per_cluster.append(time_and_cluster_per_cell["velocity_pseudotime"][time_and_cluster_per_cell["leiden"]==cluster].mean())

# Functions to Enable Downstream Analysis of Gene Sets for Single Cell Datasets

# Convert to Gene.By.Sample.Score.Matrix
def velocity_score_to_gct(adata, cluster_out, outname):
    import re
    import numpy as np
    import pandas as pd
    unique_values = set()
    for col in scv.DataFrame(adata.uns['rank_velocity_genes']['names']):
        unique_values.update(scv.DataFrame(
            adata.uns['rank_velocity_genes']['names'])[col])
    unique_values = list(unique_values)
    unique_values.sort()

    gene_by_cluster = pd.DataFrame(columns=scv.DataFrame(
        adata.uns['rank_velocity_genes']['names']).columns, index=unique_values)
    for col in scv.DataFrame(adata.uns['rank_velocity_genes']['names']):
        gene_by_cluster[col] = list(scv.DataFrame(adata.uns['rank_velocity_genes']['scores'])[
            col][np.argsort(scv.DataFrame(adata.uns['rank_velocity_genes']['names'])[col].values)])
    gene_by_cluster.index.name = "NAME"
    gene_by_cluster.index = gene_by_cluster.index.str.replace(
        '\\..*', '', regex=True)
    gene_by_cluster.insert(loc=0, column='Description', value="NA")
    text_file = open(outname + "_velocity_gene_scores_by_" +
                     cluster_out + ".gct", "w")
    text_file.write('#1.2\n')
    text_file.write(str(len(gene_by_cluster)) + "\t" +
                    str(len(gene_by_cluster.columns) - 1) + "\n")
    text_file.close()
    gene_by_cluster.to_csv(outname + "_velocity_gene_scores_by_" +
                           cluster_out + ".gct", sep="\t", mode='a')


def load_ssgsea_result(ssgsea_result):
    ssgsea_df = pd.read_csv(ssgsea_result, sep='\t', header=2, index_col=[
                            0, 1], skip_blank_lines=True)
    ssgsea_df.index = ssgsea_df.index.droplevel(1)  # Drop gene descriptions
    return ssgsea_df


# Add Clusterwise ssGSEA scores to adata.obs as a cell level score for plotting
def adata_import_ssgsea_scores(adata, ssgsea_result):
    import GeneSetAnalysisFunctions
    ssgsea_df = load_ssgsea_result(ssgsea_result)
    ssgsea_cell_df = ssgsea_df.transpose()
    ssgsea_cell_df = ssgsea_cell_df.reindex(list(adata.obs['leiden']))
    ssgsea_cell_df.index = range(len(ssgsea_cell_df.index))
    ssgsea_cell_df = ssgsea_cell_df.set_index(adata.obs.index)
    adata.obs[ssgsea_cell_df.columns] = ssgsea_cell_df[ssgsea_cell_df.columns]
    return ssgsea_cell_df


def ssgsea_plot_all(adata, ssgsea_result, basis, clusters, outname, format):  # Plotting
    import GeneSetAnalysisFunctions
    import scvelo as scv
    ssgsea_cell_df = GeneSetAnalysisFunctions.adata_import_ssgsea_scores(
        adata, ssgsea_result)
    ssgsea_sets = list(ssgsea_cell_df.columns)
    for set in ssgsea_sets:
        scv.pl.velocity_embedding_stream(adata, basis=basis, color=[
                                         set, clusters], color_map='seismic', save=set + "_" + outname + "_embedding." + format)


def ssgsea_plot_hits(adata, set_hits, ssgsea_result, basis, clusters, outname, format):  # Plotting
    import GeneSetAnalysisFunctions
    import scvelo as scv
    ssgsea_cell_df = GeneSetAnalysisFunctions.adata_import_ssgsea_scores(
        adata, ssgsea_result)
    for i in range(len(set_hits)):
        set = str(set_hits.index[i])
        scv.pl.velocity_embedding_stream(adata, basis=basis, color=[
            set, clusters], color_map='seismic', add_outline=[set_hits.iloc[i][0], set_hits.iloc[i][2]], save=set + "_Cluster_" + str(set_hits.iloc[i][0]) + "_to_" + str(set_hits.iloc[i][2]) + "_" + outname + "_embedding." + format)


# Calculate the Gene Set ES Delta pairwise for every set
def create_transition_matrix(ssgsea_result, ssgsea_cell_df, set):
    import GeneSetAnalysisFunctions
    import sys
    import pandas as pd
    ssgsea_raw_df = GeneSetAnalysisFunctions.load_ssgsea_result(ssgsea_result)
    # If we get this information from ssgsea_result, then we don't need the ssgsea_cell_df, just need the unique sets from _raw_df
    ssgsea_sets = list(ssgsea_cell_df.columns)
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
def find_candidate_transitions(adata, ssgsea_result, set, threshold, silent=False):
    import GeneSetAnalysisFunctions
    import numpy as np
    import pandas as pd
    import scvelo as scv
    import math
    # connectivities confidence
    paga_conf_df = scv.get_df(
        adata, 'paga/transitions_confidence', precision=2).T
    # paga_adj_df = scv.get_df(adata, 'paga/connectivities') # , precision=2).T # connectivities adjacency
    # paga_tree_df = scv.get_df(adata, 'paga/connectivities_tree') #, precision=2).T # connectivities subtree
    ssgsea_cell_df = GeneSetAnalysisFunctions.adata_import_ssgsea_scores(
        adata, ssgsea_result)
    set_transition = GeneSetAnalysisFunctions.create_transition_matrix(
        ssgsea_result, ssgsea_cell_df, set)
    set_transition_pass = set_transition[paga_conf_df > float(threshold)]
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
                set_transition_pass == transition_values[i])[1]) + " (Enrichment Score: " + str(float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[1]))])) + ") at PAGA transition confidence >" + str(threshold))
        set_hits.append([set, str(np.where(set_transition_pass == transition_values[i])[0]).strip("[]"), str(float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[0]))])), str(np.where(set_transition_pass == transition_values[i])[1]).strip("[]"), str(
            float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[1]))])), float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[1]))]) - float(test_set[str(int(np.where(set_transition_pass == transition_values[i])[0]))])])
    return(set_hits)


# Using the results of find_candidate_transitions keep candidates that have good directionality
def find_good_transitions(adata, ssgsea_result, threshold, silent=False):
    import GeneSetAnalysisFunctions
    ssgsea_raw_df = GeneSetAnalysisFunctions.load_ssgsea_result(ssgsea_result)
    all_sets = ssgsea_raw_df.index.to_list()
    all_set_results = []
    for set in all_sets:
        set_hits = GeneSetAnalysisFunctions.find_candidate_transitions(
            adata=adata, ssgsea_result=ssgsea_result, set=set, threshold=threshold, silent=silent)
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
    all_filtered_changes = all_positive_changes.append(all_negative_changes)
    all_filtered_changes.set_index("Gene_Set", inplace=True)
    return(all_filtered_changes)

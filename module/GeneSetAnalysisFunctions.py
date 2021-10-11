# Functions to Enable Downstream Analysis of Gene Sets for Single Cell Datasets

def velocity_score_to_gct(adata, cluster_out, outname):    # Convert to Gene.By.Sample.Score.Matrix
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
    gene_by_cluster.index = gene_by_cluster.index.str.replace('\\..*', '', regex=True)
    gene_by_cluster.insert(loc=0, column='Description', value="NA")
    text_file = open(outname + "_velocity_gene_scores_by_" +
                           cluster_out + ".gct", "w")
    text_file.write('#1.2\n')
    text_file.write(str(len(gene_by_cluster)) + "\t" + str(len(gene_by_cluster.columns) - 1) + "\n")
    text_file.close()
    gene_by_cluster.to_csv(outname + "_velocity_gene_scores_by_" +
                           cluster_out + ".gct", sep="\t", mode='a')


def load_ssgsea_scores(adata, ssgsea_result):    # Add Clusterwise ssGSEA scores to adata.obs as a cell level score for plotting
    ssgsea_df = pd.read_csv(ssgsea_result, sep='\t', header=2, index_col=[
                            0, 1], skip_blank_lines=True)
    ssgsea_df.index = ssgsea_df.index.droplevel(1)  # Drop gene descriptions
    ssgsea_df = ssgsea_df.transpose()
    ssgsea_df = ssgsea_df.reindex(list(adata.obs['leiden']))
    ssgsea_df.index = range(len(ssgsea_df.index))
    ssgsea_df = ssgsea_df.set_index(adata.obs.index)
    adata.obs[ssgsea_df.columns] = ssgsea_df[ssgsea_df.columns]
    return ssgsea_df

def ssgsea_plot(adata, ssgsea_result, basis, clusters ,outname, format):# Plotting
    import GeneSetAnalysisFunctions
    import scvelo as scv
    ssgsea_df = GeneSetAnalysisFunctions.load_ssgsea_scores(adata, ssgsea_result)
    ssgsea_sets = list(ssgsea_df.columns)
    for set in ssgsea_sets:
        scv.pl.velocity_embedding_stream(adata, basis=basis, legend_loc='right', color=[
                                         set, clusters], color_map='seismic', save=set + "_" + outname + "_embedding." + format)

def create_transition_matrix(ssgsea_result, ssgsea_df, set):
    ssgsea_sets = list(ssgsea_df.columns)
    set_transition = pd.DataFrame(columns=ssgsea_raw_df.columns, index=ssgsea_raw_df.columns)
    test_set = ssgsea_raw_df[:set]
    for first_cluster in test_set.columns:
        for second_cluster in test_set.columns:
            set_transition.at[first_cluster,second_cluster] = float(test_set[second_cluster]) - float(test_set[first_cluster])
    return set_transition

def find_outlier_transitions(adata, ssgsea_result, ssgsea_df, set):
    import GeneSetAnalysisFunctions
    set_transition = GeneSetAnalysisFunctions.create_transition_matrix(ssgsea_result, ssgsea_df, set)
    set_transition_pass = set_transition[paga_df>0]
    set_transition_pass_list = set_transition_pass.values.tolist()
    flat_set_transition_pass_list = [item for sublist in set_transition_pass_list for item in sublist if math.isnan(item) == False]
    mean = np.mean(flat_set_transition_pass_list)
    standard_deviation = np.std(flat_set_transition_pass_list)
    distance_from_mean = abs(flat_set_transition_pass_list - mean)
    max_deviations = 2
    outlier = distance_from_mean > max_deviations * standard_deviation
    outlier_loc = list(np.where(outlier)[0])
    transition_outlier_values = np.array(flat_set_transition_pass_list)[outlier_loc]
    transition_outlier_values = list(transition_outlier_values)
    for i in range(len(transition_outlier_values)):
        print(np.where(set_transition_pass==transition_outlier_values[i]))

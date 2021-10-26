
# Reimplementation of the R ssGSEA GMT Parser
# Reads a gene set database file (in GMX file format)
# and creates an Pandas Datafrme with each row corresponding to a single
# gene set and containing a list of the gene names making up
# that gene set.  Gene sets that do not satisfy the min and max threshold
# criteria will be filtered out. Returned in a dict with other information
def read_genesets_gmt(gs_db, thres_min=2, thres_max=2000):
    import pandas as pd
    import numpy as np
    with open(gs_db) as f:
        temp = f.read().splitlines()
    max_Ng = len(temp)
    # temp_size_G will contain size of each gene set
    temp_size_G = list(range(max_Ng))
    for i in range(max_Ng):
        temp_size_G[i] = len(temp[i].split("\t")) - 2
    max_size_G = max(temp_size_G)
    gs = pd.DataFrame(np.nan, index=range(max_Ng), columns=range(max_size_G))
    temp_names = list(range(max_Ng))
    temp_desc = list(range(max_Ng))
    gs_count = 0
    for i in range(max_Ng):
        gene_set_size = len(temp[i].split("\t")) - 2
        gs_line = temp[i].split("\t")
        gene_set_name = gs_line[0]
        gene_set_desc = gs_line[1]
        gene_set_tags = list(range(gene_set_size))
        for j in range(gene_set_size):
            gene_set_tags[j] = gs_line[j + 2]
        if np.logical_and(gene_set_size >= thres_min, gene_set_size <= thres_max):
            temp_size_G[gs_count] = gene_set_size
            gs.iloc[gs_count] = gene_set_tags + \
                list(np.full((max_size_G - temp_size_G[gs_count]), np.nan))
            temp_names[gs_count] = gene_set_name
            temp_desc[gs_count] = gene_set_desc
            gs_count = gs_count + 1
    Ng = gs_count
    gs_names = list(range(Ng))
    gs_desc = list(range(Ng))
    size_G = list(range(Ng))
    gs_names = temp_names[0:Ng]
    gs_desc = temp_desc[0:Ng]
    size_G = temp_size_G[0:Ng]
    gs.dropna(how='all', inplace=True)
    gs.index = gs_names
    return({'N_gs': Ng, 'gs': gs, 'gs_names': gs_names, 'gs_desc': gs_desc, 'size_G': size_G, 'max_N_gs': max_Ng})


# Reimplementation of the R ssGSEA GMX Parser
# Reads a gene set database file (in GMX file format)
# and creates an Pandas Datafrme with each row corresponding to a single
# gene set and containing a list of the gene names making up
# that gene set.  Gene sets that do not satisfy the min and max threshold
# criteria will be filtered out. Returned in a dict with other information
def read_genesets_gmx(gs_gmx, thres_min=2, thres_max=2000):
    import pandas as pd
    import numpy as np
    df_temp = pd.read_csv(
        gs_gmx, sep='\t', skip_blank_lines=True).transpose().dropna(how='all')
    all_gs_names = df_temp.index.tolist().copy()
    all_gs_desc = df_temp[0].tolist().copy()
    all_gs = df_temp.drop(labels=0, axis=1)
    all_gs_sizes = all_gs.count(axis=1).tolist()
    pass_thresholds = np.logical_and(all_gs.count(
        axis=1) >= thres_min, all_gs.count(axis=1) <= thres_max).to_list()
    gs_names = np.array(all_gs_names)[np.array(
        pass_thresholds)].tolist().copy()
    gs_desc = np.array(all_gs_desc)[np.array(pass_thresholds)].tolist().copy()
    gs_sizes = np.array(all_gs_sizes)[np.array(
        pass_thresholds)].tolist().copy()
    gs = all_gs[pass_thresholds]
    max_Ng = len(all_gs_names)
    Ng = len(gs_names)
    gs.columns = range(len(gs.columns))
    # N_gs = number of gene sets defined in gmx file that satisfy the min and max thresholds
    # gs = matrix containing gene set collections, one per line, satisfying min/max thresholds
    # gs_names = vector of names of gene sets (of length N_gs)
    # gs_desc = vector of descriptions of gene sets (of length N_gs)
    # size_G = vector with sizes of each gene set (of length N_gs)
    # max_N_gs = total number of gene sets defined in gmx file; includes those that do not satisfy min/max thresholds
    return({'N_gs': Ng, 'gs': gs, 'gs_names': gs_names, 'gs_desc': gs_desc, 'size_G': gs_sizes, 'max_N_gs': max_Ng})


def read_chip(chip):
    import os
    import sys
    import pandas as pd
    chip_df = pd.read_csv(chip, sep='\t', index_col=0, skip_blank_lines=True)
    return(chip_df)


def collapse_dataset(dataset, chip, mode="sum"):
    import GeneSetAnalysisFunctions
    import pandas as pd
    if isinstance(chip, pd.DataFrame) == False:
        chip = GeneSetAnalysisFunctions.read_chip(chip)
    joined_df = chip.join(dataset, how='inner')
    joined_df.reset_index(drop=True, inplace=True)
    annotations = joined_df[["Gene Symbol",
                             "Gene Title"]].drop_duplicates().copy()
    joined_df.drop("Gene Title", axis=1, inplace=True)
    if mode == "sum":
        collapsed_df = joined_df.groupby(["Gene Symbol"]).sum()
    if mode == "mean":
        collapsed_df = joined_df.groupby(["Gene Symbol"]).mean()
    if mode == "median":
        collapsed_df = joined_df.groupby(["Gene Symbol"]).median()
    if mode == "max":
        collapsed_df = joined_df.groupby(["Gene Symbol"]).max()
    collapsed_df.index.name = "NAME"
    return(collapsed_df)

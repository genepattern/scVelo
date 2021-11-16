# Copyright (c) 2012-2020 Broad Institute, Inc., Massachusetts Institute of Technology, and Regents of the University of California.  All rights reserved.

# ssGSEA

# processes the cmd line for the ssGSEA.project.dataset
import sys, argparse

def main():
    usage = "%prog [options]" + "\n"
    ap = argparse.ArgumentParser()
    ap.add_argument("-l", "--libdir", action="store",
                    dest="libdir", help="Working directory to load ssSGEAlib from.")
    ap.add_argument("-i", "--inputgct", action="store",
                    dest="input_gct_filename", help="Input GCT file containing expression data.")
    ap.add_argument("-o", "--outprefix", action="store",
                    dest="output_prefix", nargs='?', default=None, help="Output filename prefix.  If unspecified, is set to <prefix of input_gct_file>.PROJ.")
    ap.add_argument("-D", "--gsdb", action="store",
                    dest="gene_sets_db_list_filename", help="Gene sets database file in GMT of GMX format.")
    ap.add_argument("-c", "--column", action="store",
                    dest="gene_symbol_column", nargs='?', default='Name', help="Specifies the gct data column containing gene symbols.  In most cases this will be 'Name'.")
    ap.add_argument("-s", "--set", action="store",
                    dest="gene_set_selection", nargs='?', default='ALL', help="Comma-separated list of gene set names on which to project expression data. Set to 'ALL' if projecting to all defined gene sets.")
    ap.add_argument("-n", "--normalize", action="store",
                    dest="sample_normalization_method", nargs='?', default='none', help="Method used to normalize gene expression data. Valid options are 'none' (default), 'rank', 'log', and 'log.rank'")
    ap.add_argument("-w", "--weight", action="store",
                    dest="weighting_exponent", nargs='?', type=float, default=0.75, help="Exponent for weighting enrichment calculation. Module authors strongly recommend against changing from default.")
    ap.add_argument("-r", "--velocity_weight", action="store",
                    dest="velocity_weight", nargs='?', default=False, help="Apply a weighting factor to the data based on scores from scVelo's rank_velocity_genes")
    ap.add_argument("-v", "--min", action="store",
                    dest="min_overlap", nargs='?', type=int, default=1, help="Gene sets with overlap smaller than this are excluded from the analysis.")
    ap.add_argument("-C", "--combine", action="store",
                    dest="combine_mode", nargs='?', default='combine.add', help="Options for combining enrichment scores for paired *_UP and *_DN gene sets. Valid options are 'combine.off' (default, do not combine paired sets), 'combine.replace' (replace the parent gene sets with just the combined set). 'combine.add' (report the combined set in addition to the parent sets).")
    ap.add_argument("-m", "--collapse_method", action="store",
                    dest="collapse_method", nargs='?', default='none', help="Method to use for collapsing dataset from gene ids to gene symbols. Valid options are 'none' (default, do not collapse the dataset), 'sum' (recommended for RNA-seq data), 'max' (GSEA Desktop default, recommended for microarrays), 'mean', and 'median'.")
    ap.add_argument("-f", "--chip_file", action="store",
                    dest="chip_file", nargs='?', default=None, help="CHIP file containing Gene ID to Gene Symbol mappings for use with collapse_method!='none'.")
    ap.add_argument("-t", "--genes_min_nonzero_cells_threshold", action="store",
                    dest="cell_threshold", nargs='?', default=0, help="Sets a threshold for how many cells must express a gene for it to be included in the output matrix.")

    options = ap.parse_args()

    sys.path.insert(1, options.libdir)
    import ssGSEAlib

    if options.input_gct_filename == None:
        sys.exit("Missing input_gct_filename")

    if options.input_gct_filename == ssGSEAlib.check_extension(options.input_gct_filename,".h5ad"):
        import GeneSetAnalysisFunctions
        import scanpy as sc
        adata = sc.read(options.input_gct_filename)
        if options.output_prefix == None:
            temp = options.input_gct_filename.split("/") # Extract input file name
            input_file_name = temp[-1]
            temp = input_file_name.split(".h5ad")
        else:
            output_ds = options.output_prefix + ".gct"
        if options.gene_symbol_column == "pseudobulk_counts":
            options.input_gct_filename = GeneSetAnalysisFunctions.make_pseudobulk(adata,genes_min_nonzero_cells = int(options.cell_threshold), outname = temp[0], velocity_weight=GeneSetAnalysisFunctions.str_to_bool(options.velocity_weight))
            temp[0] = options.input_gct_filename["outname"]
        elif options.gene_symbol_column.upper() == "NAME":
            print("Please use one of the single-cell specific options. Falling back to pseudobulk count data.")
            options.input_gct_filename = GeneSetAnalysisFunctions.make_pseudobulk(adata,genes_min_nonzero_cells = int(options.cell_threshold), outname = temp[0], velocity_weight=GeneSetAnalysisFunctions.str_to_bool(options.velocity_weight))
            temp[0] = options.input_gct_filename["outname"]
        else:
            options.input_gct_filename = GeneSetAnalysisFunctions.get_gene_values(adata,key=options.gene_symbol_column, genes_min_nonzero_cells = int(options.cell_threshold), outname = temp[0], velocity_weight=options.velocity_weight)
            temp[0] = options.input_gct_filename["outname"]
        output_ds =  temp[0]+".PROJ.gct"
        options.gene_symbol_column = "Name"
    else:
        if options.output_prefix == None:
            temp = options.input_gct_filename.split("/") # Extract input file name
            input_file_name = temp[-1]
            temp = input_file_name.split(".gct")
            output_ds =  temp[0]+".PROJ.gct"
        else:
            output_ds = options.output_prefix + ".gct"

    if options.gene_sets_db_list_filename != None:
        with open(options.gene_sets_db_list_filename) as f:
            gene_sets_dbfile_list = f.read().splitlines()
    else:
        sys.exit("No Gene Set DB files provided")

    if options.collapse_method != "none":
        input_ds = ssGSEAlib.collapse_dataset(options.input_gct_filename, options.chip_file, method=options.collapse_method)
    else:
        input_ds = options.input_gct_filename

    ssGSEAlib.ssGSEA_project_dataset(input_ds,
                                            output_ds,
                                            gene_sets_dbfile_list,
                                            gene_symbol_column = options.gene_symbol_column,
                                            gene_set_selection = options.gene_set_selection,
                                            sample_norm_type = options.sample_normalization_method,
                                            weight = options.weighting_exponent,
                                            combine_mode = options.combine_mode,
                                            min_overlap = options.min_overlap)

if __name__ == '__main__':
    main()

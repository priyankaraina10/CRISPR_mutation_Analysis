import itertools
import argparse
import logging
import sys
import numpy as np
import pandas as pd
import scipy.stats

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def parse_commandline():
    """Parse command line parameters.
    Return type: Object
    """
    parser = argparse.ArgumentParser(description="Association of models and mutations")
    parser.add_argument("--mutation_data", metavar="FILE.tsv", help="mutation data file", required=True)
    parser.add_argument("--gene_data", metavar="FILE.tsv", help="gene data file", required=True)
    parser.add_argument("--all_mutations", metavar="FILE.txt", help="all mutations data file", required=True)
    parser.add_argument("--all_genes", metavar="FILE.txt", help="all_genes data file", required=True)
    parser.add_argument("--models_mutations", metavar="FILE.txt", help="models mutations data file", required=True)
    parser.add_argument("--models_genes", metavar="FILE.txt", help="models genes data file", required=True)
    parser.add_argument("--models_common", metavar="FILE.txt", help="models common data file", required=True)
    parser.add_argument("--end_mutation_data", metavar="FILE.tsv", help="end mutation data file", required=True)
    parser.add_argument("--end_gene_data", metavar="FILE.tsv", help="end gene data file", required=True)
    parser.add_argument("--output", metavar="FILE.tsv", help="Output file name", required=True)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()

def main():
    params = parse_commandline()

    #load gene data
    logging.info("Reading gene and mutation data files")
    try:
        mutation_data = pd.read_csv(params.mutation_data, sep="\t", index_col="Mutation")
        gene_data = pd.read_csv(params.gene_data, sep="\t", index_col="Model")
    except pd.errors.ParserError as pe:
        logging.info("Can't read input data files with the right format. Creating empty default file...")
        sys.exit(0)


    #load mutations and genes present in gene and mutation data files
    logging.info("Reading mutations and genes present in gene and mutation data files")
    all_mutations = [line.rstrip() for line in open(params.all_mutations)]
    all_genes = [line.rstrip() for line in open(params.all_genes)]

    #Models in genedata and mutation_data file
    logging.info("Reading models in gene and mutation data files")
    models_mutations = [line.rstrip() for line in open(params.models_mutations)]
    models_genes = [line.rstrip() for line in open(params.models_genes)]

    ##models present in both files
    models_common = [line.rstrip() for line in open(params.models_common)]

    ### To know whether the presence of that variant is associated to a change in the outcome of the CRISPR knock-out experiment for that gene, across all cell lines.
    ## get only the common models because those are the only ones we can compare
    end_mutation_table = pd.read_csv(params.end_mutation_data,
                 sep = "\t", index_col ="Models")

    end_gene_table = pd.read_csv(params.end_gene_data,
                 sep = "\t", index_col ="Models")

    #calculate p_values
    all_p_values = []
    genes = []
    mutations = []
    for gene, mutation in itertools.product(all_genes, all_mutations):
        gene_table_to_test = np.array(end_gene_table[gene])
        mutation_table_to_test = np.array(end_mutation_table[mutation])
        p_value = scipy.stats.wilcoxon(gene_table_to_test, mutation_table_to_test)[1]
        genes.append(gene)
        mutations.append(mutation)
        all_p_values.append(p_value)

    ## calculating adjusted_pvalues
    adjusted_pvalues = p_adjust_bh(all_p_values)
    d = {'Gene': genes, 'Mutation': mutations, 'pvalue': all_p_values,'p_adj':adjusted_pvalues}
    df = pd.DataFrame(data=d)

    #Print output file
    logging.info("Writing into file")
    df.to_csv(params.output, sep="\t")

if __name__ == "__main__":
    main()



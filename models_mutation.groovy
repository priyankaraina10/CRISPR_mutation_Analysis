#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Defining data directories
 */

params.data_dir = "home/ubuntu/CRISPR_Mutation_Analysis/"
params.out_dir = "/home/ubuntu/CRISPR_Mutation_Analysis/results/"

workflow {
    association_mutation_models(params.data_dir)
}

workflow association_mutation_models {
    take:
        data_dir

    main:
        // Input data
        mutation_data_ch = Channel.fromPath(params.data_dir + '/test/Mutations.tsv').first()
        gene_data_ch = Channel.fromPath(params.data_dir + '/test/Gene_KOs.tsv').first()
        all_mutations_ch = Channel.fromPath(params.data_dir + '/test/all_mutations.txt').first()
        all_genes_ch = Channel.fromPath(params.data_dir + '/test/all_genes.txt').first()
        models_mutations_ch = Channel.fromPath(params.data_dir + '/test/models_mutations.txt').first()
        models_genes_ch = Channel.fromPath(params.data_dir + '/test/models_genes.txt').first()
        models_common_ch = Channel.fromPath(params.data_dir + '/test/models_common.txt').first()
        end_mutation_data_ch = Channel.fromPath(params.data_dir + '/test/end_mutation_table.tsv').first()
        end_gene_data_ch = Channel.fromPath(params.data_dir + '/test/end_gene_table.tsv').first()
        association_mm(mutation_data_ch, gene_data_ch, all_mutations_ch, all_genes_ch, models_mutations_ch, models_genes_ch, models_common_ch, end_mutation_data_ch, end_gene_data_ch)
}

/*
 * association_mutation_models analysis
 */

process association_mm {
    container 'CRISPR_Mutation_Analysis'
    storeDir params.out_dir  + 'association_mm'
    input:
        path("Mutations.tsv")
        path("Gene_KOs.tsv")
        path("all_mutations.txt")
        path("all_genes.txt")
        path("models_mutations.txt")
        path("models_genes.txt")
        path("models_common.txt")
        path("end_mutation_table.tsv")
        path("end_gene_table.tsv")

    output:
        path("association_mm.tsv")

    script:
        """
           python3 /association_mutation_models.py --mutation_data Mutations.tsv --gene_data Gene_KOs.tsv --all_mutations all_mutations.txt --all_genes all_genes.txt --models_mutations models_mutations.txt --models_genes models_genes.txt --models_common models_common.txt --end_mutation_data end_mutation_table.tsv --end_gene_data end_gene_table.tsv --output association_mm.tsv
        """
}

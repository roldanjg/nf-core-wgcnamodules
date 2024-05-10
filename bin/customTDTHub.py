#!/usr/bin/env python3

# Written by Joaquin Grau Roldán . Released under the MIT license.

import argparse
import sys
import pandas as pd

def parse_args(args=None):
    Description = "Reformat nf-core/rnaseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--genes_modules", type=str, required=True, help="")
    return parser.parse_args(args)

def load_WGCNA_gene_modules(wgcnaOut=None):

    gene_info_df = pd.read_csv(wgcnaOut, sep=',')
    print()
    gene_info_df = gene_info_df.set_index('geneid')
    pval_color = []
    prob_color = []
    gene_info_df = gene_info_df[gene_info_df.moduleColor != 0]
    for _, row in gene_info_df.iterrows():
        pval_color.append(row['p.MM.'+(str(int(row.moduleColor)))])
        prob_color.append(row['MM.'+(str(int(row.moduleColor)))])
    gene_info_df['moduleColor_prob'] = prob_color
    gene_info_df['moduleColor_pval'] = pval_color
    ontology_info_df = gene_info_df[['moduleColor','moduleColor_prob','moduleColor_pval',]]
    gene_info_df = gene_info_df.drop(columns=['moduleColor','moduleColor_pval','moduleColor_prob'])
    pvals = []
    probs = []
    for column in gene_info_df.columns:
        if 'p.' in column:
            pvals.append(column)
        if not 'p.' in column and 'MM' in column:
            probs.append(column)

    ontology_info_df['max_prob_module'] = gene_info_df[probs].idxmax(axis=1)
    ontology_info_df['max_prob_value'] = gene_info_df[probs].max(axis=1)
    ontology_info_df['min_pval_module'] = gene_info_df[pvals].idxmin(axis=1)
    ontology_info_df['min_pval_value'] = gene_info_df[pvals].min(axis=1)
    return ontology_info_df

def generate_individual_module_genes(df):

    for module in df['max_prob_module'].unique():
        colorgenes = df.loc[
            df['max_prob_module'] == module
        ]
        with open(f'{module.split(".")[1]}.csv','w') as save_f: #module
            for line_gene in list(colorgenes.index):
                save_f.write(f'{line_gene}\n')

def main(args=None):
    args = parse_args(args)
    modules_df = load_WGCNA_gene_modules(args.genes_modules)
    generate_individual_module_genes(modules_df)
if __name__ == "__main__":
    sys.exit(main())

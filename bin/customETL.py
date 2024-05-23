#!/usr/bin/env python3

# Written by Joaquin Grau Roldán . Released under the MIT license.
import argparse
from math import log2
import pandas as pd
import sys

def parse_args(args=None):
    Description = ""
    Epilog = ""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--input_wgcna", type=str, required=True, help="")
    parser.add_argument("--contrast_wgcna", type=str, required=False, help="")
    parser.add_argument("--tpms_wgcna", type=str, required=False, help="")
    parser.add_argument("--norm", default="tpm", type=str, required=False, help="")
    parser.add_argument("--genes", default="no", type=str, required=False, help="")
    parser.add_argument("--fdr", type=float, help="")
    parser.add_argument("--log2ratio", type=float, help="")

    return parser.parse_args(args)

def genes_extraction(genes, fdr,log2fold):
    print(fdr,log2fold)
    path_list = genes.split('+c+')
    genes_list = set()
    for genes_path in path_list:
        with open(genes_path, 'r') as genes_to_filter:
            next(genes_to_filter)
            for line in genes_to_filter:
                gene_line = line.strip().split('\t')
                gene_id = gene_line[0]
                gene_log2fold = gene_line[2]
                gene_fdr = gene_line[5]
                if gene_fdr == 'NA':
                    continue
                if (
                    abs(
                        float(gene_log2fold)
                        ) > log2fold and float(gene_fdr) < fdr
                ):
                    genes_list.add(gene_id)

    with open("diff_selected_genes.txt", "w") as genes_to_save:
        for gene in genes_list:
            genes_to_save.write(f"{gene}\n")

def load_data(input_wgcna,contrast_wgcna,tpms_wgcna,genes):
    metadata_dict = {'samples':{},'comparations':{}, 'genes':[]}
    with open(input_wgcna, "r", encoding="utf-8-sig") as fin:
        next(fin)
        for line in fin:
            line_sp = line.split(',')
            category = metadata_dict['samples'].setdefault(line_sp[1], [])
            category.append(line_sp[0])
    with open(contrast_wgcna, "r", encoding="utf-8-sig") as fin:
        next(fin)
        for line in fin:
            line_sp = line.strip().split(',')
            metadata_dict['comparations'][line_sp[0]] = {'reference':line_sp[2],'target':line_sp[3]}
    if genes.endswith('.txt'):
        with open(genes, "r", encoding="utf-8-sig") as fin:
            for line in fin:
                metadata_dict['genes'].append(line.strip())

    tpms = pd.read_csv(tpms_wgcna, sep='\t').drop(columns='gene_name')
    return metadata_dict, tpms

def calculate_mean(metadata, tpmsdf):
    tpmsdf = tpmsdf.set_index("gene_id")
    tpms_cols = []
    for metadata_condition in metadata['samples']:
        replicates = metadata['samples'][metadata_condition]
        tpms_cols.append(metadata_condition)
        tpmsdf[metadata_condition] = tpmsdf[replicates].mean(axis=1)
    return tpmsdf[tpms_cols]

def select_specific_genes(metadata, tpmsdf):
    return tpmsdf[tpmsdf['gene_id'].isin(metadata['genes'])]

def calculate_ratio(metadata, tpmsdf):

    ratio_cols = []
    tpmsdf = tpmsdf.map(lambda c: log2(c + 1))
    for comparation in metadata['comparations']:
        ratio_cols.append(comparation)
        target = metadata['comparations'][comparation]['target']
        reference = metadata['comparations'][comparation]['reference']
        tpmsdf[comparation] = tpmsdf[target] - tpmsdf[reference]

    return tpmsdf[ratio_cols]

def main(args=None):
    args = parse_args(args)
    # generate genes.txt from diff expression
    if 'tsv+c+' in args.genes:
        genes_extraction(args.genes, args.log2ratio, args.fdr)
        args.genes = 'diff_selected_genes.txt'

    # extract
    metadata, tpmsdf = load_data(
        args.input_wgcna,args.contrast_wgcna,args.tpms_wgcna,args.genes
        )
    # transform

    if args.genes.endswith('.txt'):
        tpmsdf = select_specific_genes(metadata, tpmsdf)

    tpmsdf = calculate_mean(metadata, tpmsdf)
    if args.norm == 'ratio':
        tpmsdf = calculate_ratio(metadata, tpmsdf)
    # load
    tpmsdf.reset_index().T.to_csv("WGCNA_input_filtered.tpms.csv", sep=",", header=False)


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3

# Written by Joaquin Grau Roldán . Released under the MIT license.
import argparse
from math import log2
import pandas as pd
import sys

def parse_args(args=None):
    Description = "Reformat nf-core/rnaseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--input_wgcna", type=str, required=True, help="")
    parser.add_argument("--contrast_wgcna", type=str, required=False, help="")
    parser.add_argument("--tpms_wgcna", type=str, required=False, help="")
    parser.add_argument("--norm", default="tpm", type=str, required=False, help="")
    parser.add_argument("--genes", default="no", type=str, required=False, help="")
    parser.add_argument("--fdr", type=float, help="")
    parser.add_argument("--log2ratio", type=float, help="")

    return parser.parse_args(args)

# def extract():

# def make_dir(path):
#     if len(path) > 0:
#         try:
#             os.makedirs(path)
#         except OSError as exception:
#             if exception.errno != errno.EEXIST:
#                 raise exception


# def print_error(error, context="Line", context_str=""):
#     error_str = f"ERROR: Please check samplesheet -> {error}"
#     if context != "" and context_str != "":
#         error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
#     print(error_str)
#     sys.exit(1)


# def check_samplesheet(file_in, file_out):
#     """
#     This function checks that the samplesheet follows the following structure:

#     sample_name,replicate,control_id,class,stress_type
#     SAMPLE_PE,1,control,heat
#     SAMPLE_PE,2,control,heat
#     SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,,forward

#     For an example see:
#     https://github.com/nf-core/test-datasets/blob/rnaseq/samplesheet/v3.1/samplesheet_test.csv
#     """

#     sample_mapping_dict = {}
#     with open(file_in, "r", encoding="utf-8-sig") as fin:

#         ## Check header
#         MIN_COLS = 3
#         HEADER = ["sample", "fastq_1", "fastq_2", "strandedness"]
#         header = [x.strip('"') for x in fin.readline().strip().split(",")]
#         if header[: len(HEADER)] != HEADER:
#             print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
#             sys.exit(1)

#         ## Check sample entries
#         for line in fin:
#             if line.strip():
#                 lspl = [x.strip().strip('"') for x in line.strip().split(",")]

#                 ## Check valid number of columns per row
#                 if len(lspl) < len(HEADER):
#                     print_error(
#                         f"Invalid number of columns (minimum = {len(HEADER)})!",
#                         "Line",
#                         line,
#                     )

#                 num_cols = len([x for x in lspl[: len(HEADER)] if x])
#                 if num_cols < MIN_COLS:
#                     print_error(
#                         f"Invalid number of populated columns (minimum = {MIN_COLS})!",
#                         "Line",
#                         line,
#                     )

#                 ## Check sample name entries
#                 sample, fastq_1, fastq_2, strandedness = lspl[: len(HEADER)]
#                 if sample.find(" ") != -1:
#                     print(f"WARNING: Spaces have been replaced by underscores for sample: {sample}")
#                     sample = sample.replace(" ", "_")
#                 if not sample:
#                     print_error("Sample entry has not been specified!", "Line", line)

#                 ## Check FastQ file extension
#                 for fastq in [fastq_1, fastq_2]:
#                     if fastq:
#                         if fastq.find(" ") != -1:
#                             print_error("FastQ file contains spaces!", "Line", line)
#                         if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz") and not fastq.endswith(".fastq"):
#                             print_error(
#                                 "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
#                                 "Line",
#                                 line,
#                             )

#                 ## Check strandedness
#                 strandednesses = ["unstranded", "forward", "reverse", "auto"]
#                 if strandedness:
#                     if strandedness not in strandednesses:
#                         print_error(
#                             f"Strandedness must be one of '{', '.join(strandednesses)}'!",
#                             "Line",
#                             line,
#                         )
#                 else:
#                     print_error(
#                         f"Strandedness has not been specified! Must be one of {', '.join(strandednesses)}.",
#                         "Line",
#                         line,
#                     )

#                 ## Auto-detect paired-end/single-end
#                 sample_info = []  ## [single_end, fastq_1, fastq_2, strandedness]
#                 if sample and fastq_1 and fastq_2:  ## Paired-end short reads
#                     sample_info = ["0", fastq_1, fastq_2, strandedness]
#                 elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
#                     sample_info = ["1", fastq_1, fastq_2, strandedness]
#                 else:
#                     print_error("Invalid combination of columns provided!", "Line", line)

#                 ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, strandedness ]]}
#                 sample_info = sample_info + lspl[len(HEADER) :]
#                 if sample not in sample_mapping_dict:
#                     sample_mapping_dict[sample] = [sample_info]
#                 else:
#                     if sample_info in sample_mapping_dict[sample]:
#                         print_error("Samplesheet contains duplicate rows!", "Line", line)
#                     else:
#                         sample_mapping_dict[sample].append(sample_info)

#     ## Write validated samplesheet with appropriate columns
#     if len(sample_mapping_dict) > 0:
#         out_dir = os.path.dirname(file_out)
#         make_dir(out_dir)
#         with open(file_out, "w") as fout:
#             fout.write(
#                 ",".join(["sample", "single_end", "fastq_1", "fastq_2", "strandedness"] + header[len(HEADER) :]) + "\n"
#             )
#             for sample in sorted(sample_mapping_dict.keys()):

#                 ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
#                 if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
#                     print_error(
#                         f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
#                         "Sample",
#                         sample,
#                     )

#                 ## Check that multiple runs of the same sample are of the same strandedness
#                 if not all(x[3] == sample_mapping_dict[sample][0][3] for x in sample_mapping_dict[sample]):
#                     print_error(
#                         f"Multiple runs of a sample must have the same strandedness!",
#                         "Sample",
#                         sample,
#                     )

#                 for idx, val in enumerate(sample_mapping_dict[sample]):
#                     fout.write(",".join([f"{sample}_T{idx+1}"] + val) + "\n")
#     else:
#         print_error(f"No entries to process!", "Samplesheet: {file_in}")

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
    print(metadata['comparations'])
    ratio_cols = []
    tpmsdf = tpmsdf.applymap(lambda c: log2(c + 1))
    for comparation in metadata['comparations']:
        ratio_cols.append(comparation)
        target = metadata['comparations'][comparation]['target']
        reference = metadata['comparations'][comparation]['reference']
        tpmsdf[comparation] = tpmsdf[target] - tpmsdf[reference]

    return tpmsdf[ratio_cols]

def main(args=None):
    args = parse_args(args)
    # extract
    metadata, tpmsdf = load_data(args.input_wgcna,args.contrast_wgcna,args.tpms_wgcna,args.genes)
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
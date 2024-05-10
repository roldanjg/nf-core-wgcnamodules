nextflow run  nf-core-wgcnamodules --input <path>/samplesheet_wgcna.csv --contrast <path>contrasts_wgcna.csv --outdir <path> --salmon_dir <path> -profile conda --norm "tpm" --diffgenes false --fdr 0.05 --log2ratio 1  --max_cpus 6 --max_memory 7.GB
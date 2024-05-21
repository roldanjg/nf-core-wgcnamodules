# nf-core-wgcnamodules: Usage

## Introduction

The computational pipeline detailed in this documentation one of the steps in a 3 steps workflow to  to infer TF regulators applicable to 60 plant species from RNA sequencing (RNA-seq) data as the starting point.It consists of 3 main parts. The first part, is RNA-Seq quantification, using [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) pipeline. The second part, nf-core-wgcnamodules, is executed from [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) results and consists of an optional idetification of differentially expressed genes between pairs of contrasts, and the use of gene expression levels to perform a clustering analysis using Weighted Correlation Network Analysis (WGCNA). In the last step, genes corresponding to each cluster are used for the identification of enriched TF binding sites that could help to determine the most relevant transcriptional regulators involved in the biological process under study using the webpage [`TDTHub`](http://acrab.cnb.csic.es/TDTHub/). This documentation detail the use of nf-core-wgcnamodules.
<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet inputs

You will need to create two samplesheets with information about the samples you would like to analyse before running the pipeline, and the contrasts you would like to study between those samples. Both have to be a comma-separated file with the number of columns  and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]' --contrast '[path to samplesheet file]'
```

### Input samplesheet

You can adapt the shamplesheet input file from [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) as the pipelines are intended to run in a complementary way. It is important that the `sample` field matches the header in the Salmon folder files that are generated from [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq).


```csv title="samplesheet.csv"
sample,condition,replicate
CONTROL1_REP1,CONTROL1,1
CONTROL1_REP2,CONTROL1,2
TREATMENT1_REP1,TREATMENT1,1
TREATMENT1_REP2,TREATMENT1,2
TREATMENT2_REP1,TREATMENT2,1
TREATMENT2_REP2,TREATMENT2,2
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. Same name as the 'samplesheet_rnaseq.csv' sample name from [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq). |
| `condition` | Name of the treatment, genotype or group that defines an experimental condition with one or multiple replicates.                                                             |
| `replicate` |Number of the biological replicate. |

An [example input samplesheet](../assets/samplesheet_wgcna.csv) has been provided with the pipeline.

### Contrast samplesheet

The contrasts file references the observations file to define groups of samples to compare. For example, based on the sample sheet above we could define contrasts like:

```csv title="contrast_wgcna.csv"
contrast,variable,control,target
TREATMENT1_vs_CONTROL1,condition,CONTROL1,TREATMENT1
TREATMENT2_vs_CONTROL1,condition,CONTROL1,TREATMENT2 
```
| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `contrast`  | A custom name used to identify the contrast. |
| `variable` | The name of the column from 'samplesheet_wgcna.csv' file that contains the condition ids. |
| `control` | The base/reference level for the contrast.  |
| `target` | The target/ non-reference level for the comparison.  

An [example contrast samplesheet](../assets/contrasts_wgcna.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core-wgcnamodules \
    -profile conda \
    --input samplesheet_wgcna.csv \ 
    --contrast contrasts_wgcna.csv \
    --salmon_dir <PATH_TO_NF-CORE/RNASEQ_SALMON_FOLDER>/salmon \
    --diff_exp_genes true \
    --outdir <OUTDIR>
```

This will launch the pipeline with the `conda` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:



```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core-wgcnamodules -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet_wgcna.csv'
contrast: './contrasts_wgcna.csv'
salmon_dir: '<PATH_TO_NF-CORE/RNASEQ_SALMON_FOLDER>/salmon' 
diff_exp_genes: true
outdir: <OUTDIR>
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Main arguments
### `--salmon_dir`
Use this to specify the path tho the location of your [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) folder with Salmon output, including raw quantified reads and TPM- (tags per million)-normalized expression matrix:

```bash
--salmon_dir '<PATH_TO_NF-CORE/RNASEQ_SALMON_FOLDER>/salmon'
```

Please note the following requirements:

1. The format of the files inside the folder should be the same as the format obtained by running salmon with [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) pipeline.
2. The folder should contain at leats these three files inside: `salmon.merged.gene_counts.tsv`, `salmon.merged.gene_lengths.tsv`, `salmon.merged.gene_tpm.tsv`.



### `--diff_exp_genes [true/false]`
Whether to perform differential expression analysis `true` or not `false` to select only DE genes for clustering. :

```bash
--diff_exp_genes true
```

It is not compatible with the option `--genes`.

### `--fdr [number] --log2fold [number]`
If `--diff_exp_genes true`, you can also select the Log2Fold and FDR threshold to extract DE genes for clustering. Default values are:

```bash
--fdr 0.05 --log2fold 1
```

### `--genes *.txt`
If `--diff_exp_genes false`, you can introduce your custom list of genes to perform WGCNA clustering.:

```bash
--genes <FILE>.txt
```
Please note the following requirements:

1. It is important to introduce a minimum of 50 genes to have at least 2 clusters. The process may fail because there are too few genes to group together.
2. The genes should appear in a `.txt` file with one gene on each line.
3. The gene names must be the same as those that appear in the `gene_id` column of the `salmon.merged.gene_tpm.tsv` file.
4. If you do not provide the `--genes *.txt` and set the `--diff_exp_genes false`, the clustering will be done with all the genes in the species genome.

An [example genes file](../assets/genes.txt) has been provided with the pipeline.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```


<!-- ### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/wgcnamodules
``` -->

<!-- ### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/wgcnamodules releases page](https://github.com/nf-core/wgcnamodules/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
::: -->

## Core Nextflow arguments (for advance users)

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,conda` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

At this time the pipeline only supports one profile: 

<!--TODO GENERATE TEST PROFILE - `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters -->

- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration (for advance users)

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).



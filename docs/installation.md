# nf-core-wgcnamodules Installation

To start using the nf-core-wgcnamodule pipeline, there are 2 steps described below:

1. [Install Nextflow and nf-core](#install-nextflow)
2. [Install the pipeline](#install-the-pipeline)

## 1) Install NextFlow

If you have conda, you can install nf-core and nexflow with just one command:

```bash
conda create --name nf-core python=3.11 nf-core=2.13.1  nextflow=23.10.1
```

**You need NextFlow version >= 23.10.1 to run this pipeline.**

See [Installation of nf-core dependencies](https://nf-co.re/docs/usage/getting_started/installation) and [nf-core](https://nf-co.re/) for further instructions on how to install and configure Nextflow for nf-core.

## 2) Install the Pipeline

The pipeline itself is not in the nf-core repository for now, so you need to install it by running:

```bash
git clone https://github.com/roldanjg/nf-core-wgcnamodules.git
```

# Differential Expression Analysis of Amblyomma americanum Data

This repository contains the workflow for performing differential expression analysis on publicly available *Amblyomma americanum* data.

## Overview of the project and the results

### Background
 
The data has been sourced from various SRA studies, and processed to quantify transcripts and perform differential expression analysis.

### Methods

#### Data Sources

The data analyzed in this repository are obtained from the following SRA studies:

- **SRP091404**: 6 samples, Single-end short reads, [Reference](https://doi.org/10.1016/j.aspen.2018.05.009)
- **SRP446981**: 24 samples, Paired-end short-read [Reference](https://doi.org/10.3389/fcimb.2023.1236785)
- **SRP032795**: 12 samples (biopanned excluded), Paired-end short-read [Reference](https://doi.org/10.1186/1471-2164-15-518)
- **SRP051699**: 4 samples, Paired-end short-read [Reference](https://doi.org/10.1371/journal.pone.0131292)
- **SRP052078; SRP052091; SRP052108; SRP052106; SRP052114; SRP052123; SRP052145; SRP052154**: 8 samples, Paired-end short-read (not published, no reference)

#### Data preprocessing

#### Counting genes

We quantified transcripts using `salmon`. 
We used the reference transcriptome assembled in the repository [2023-amblyomma-americanum-txome-assembly](https://github.com/Arcadia-Science/2023-amblyomma-americanum-txome-assembly/) to quantify read counts.
Salmon produces transcript counts however differential expression results are more accurate when gene counts, not transcript counts, are compared.
The most common way to summarize transcript counts to gene counts is to use a transcript to gene mapping file.
The R packate `tximport` uses a `tx2gene` file to sum the counts for all transcripts that encode the same gene and to report the gene-level counts.  
To generate a transcript -> gene map (`tx2gene` file) for gene-level quantification, we first mapped the reference transcripts back to the genome using uLTRA.
Then, we assigned a gene name to a transcript when it overlapped part of the genes interval as annotated in the reference GFF3 annotation file produced in the repository [2023-amblyomma-annotation](https://github.com/Arcadia-Science/2023-amblyomma-annotation).
This step is performed by the script [assign_mapped_transcripts_to_gene_by_gtf_overlap.py](./scripts/assign_mapped_transcripts_to_gene_by_gtf_overlap.py).

#### Building differential expression models 

We imported transcript counts and summarized them to gene counts using the `tximport` package function `tximport()` with the parameter `type = salmon`.
We then built differential expression models using the `DESeq2` package commands `DESeqDataSetFromTximport()` and `DESeq()`
This step is performed by the scripts [build_diffex_models.R](./scripts/build_diffex_models.R).

#### Shiny application to explore the differential expression results

The previous steps are orchestrated by a [Snakefile](./Snakefile).
The main outputs of this workflow become the inputs to the [Shiny App](./shiny).
The Shiny App hosts a reactive user interface to explore the results from the differential expression models.
Please refer to the [README in the shiny app folder](./shiny/README.md) for more information.
 
### Results

#### Batch effects limit which samples we can compare

While we had 54 input samples, we were only able to analyze 20.
An initial exploration of sample similarity is available in the notebook [20231013-differential-expression-by-groups.ipynb](./notebooks/20231013-differential-expression-by-groups.ipynb).
Below we discuss exclusion criteria for each set of samples that we removed before performing differential expression analysis.

- **SRP091404**: From the six initial samples, we excluded four, samples Um, Uf, Im, and If as these samples were all exposed to the pathogen Ehrlichia chaffeensis [Reference](https://doi.org/10.1016/j.aspen.2018.05.009). Since no other RNA-seq samples in other studies were exposed to this pathogen, we had no way to account for this variable. Further, these samples didn't have replicates, which are required to perform differential expression.
- **SRP446981**: These 24 samples were excluded because of batch effects. All samples from the study clustered tightly together and away from other similar biological replicates from other studies, indicating the batch effects were too strong to make cross-study comparisons. These samples were all unfed female whole ticks injected either with either PBS (control) or *Escherichia coli*, so it is possible that the injections caused a biological impact that led to the batch effects, but this is not possible to evaluate with the available data [Reference](https://doi.org/10.3389/fcimb.2023.1236785).
- **SRP052078; SRP052091; SRP052108; SRP052106; SRP052114; SRP052123; SRP052145; SRP052154**: These eight samples originate from two *A. americanum* cell lines. The samples all cluster tightly together and away from other samples. However, since we have no other cell line data from other studies, there is no way for us to evaluate whether they cluster alone because they have different expression or strong batch effects.

The final set of samples we used to perform differential expression are described in the table below.

|library_name   |run_accession |bioproject  |sex    |tissue         |blood_meal_hour |blood_meal_hour_range | total_spots|publication_doi                                                                     |sex_tissue              |sex_tissue_blood_meal_hour       |
|:--------------|:-------------|:-----------|:------|:--------------|:---------------|:---------------------|-----------:|:-----------------------------------------------------------------------------------|:-----------------------|:--------------------------------|
|AmbamSG7-11d   |SRR1740611    |PRJNA218793 |female |salivary_gland |168;192;216;264 |168_264               |    44443612|https://doi.org/10.1371/journal.pone.0131292                                        |female_x_salivary_gland |NA                               |
|AmbamSG72-144h |SRR1740609    |PRJNA218793 |female |salivary_gland |72;120;144      |72_144                |    49791091|https://doi.org/10.1371/journal.pone.0131292                                        |female_x_salivary_gland |female_x_salivary_gland_x_72_144 |
|AmbameSG12-18h |SRR1740608    |PRJNA218793 |female |salivary_gland |12;18;24;36;48  |12_48                 |    25095072|https://doi.org/10.1371/journal.pone.0131292                                        |female_x_salivary_gland |female_x_salivary_gland_x_12_48  |
|AmbamSGunfed   |SRR1740607    |PRJNA218793 |female |salivary_gland |0               |0                     |    53124914|https://doi.org/10.1371/journal.pone.0131292                                        |female_x_salivary_gland |NA                               |
|PL17           |SRR1027761    |PRJNA226980 |male   |whole          |72              |72_144                |      413856|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |male_x_whole            |male_x_whole_x_72_144            |
|PL3            |SRR1027763    |PRJNA226980 |male   |whole          |72              |72_144                |      850276|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |male_x_whole            |male_x_whole_x_72_144            |
|PL2            |SRR1027762    |PRJNA226980 |female |whole          |96              |72_144                |     1031591|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_whole          |female_x_whole_x_72_144          |
|PL1            |SRR1027751    |PRJNA226980 |female |whole          |24              |12_48                 |      740878|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_whole          |NA                               |
|AAFM           |SRR1027485    |PRJNA226980 |male   |whole          |72              |72_144                |       50887|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |male_x_whole            |male_x_whole_x_72_144            |
|AAUM           |SRR1027483    |PRJNA226980 |male   |whole          |0               |0                     |       78212|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |male_x_whole            |NA                               |
|AAFF           |SRR1027481    |PRJNA226980 |female |whole          |96              |72_144                |       95726|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_whole          |female_x_whole_x_72_144          |
|AAUF           |SRR1027479    |PRJNA226980 |female |whole          |0               |0                     |      352313|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_whole          |NA                               |
|MG120          |SRR1027477    |PRJNA226980 |female |midgut         |120             |72_144                |       81897|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_midgut         |female_x_midgut_x_72_144         |
|MG96           |SRR1027476    |PRJNA226980 |female |midgut         |96              |72_144                |       40713|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_midgut         |female_x_midgut_x_72_144         |
|MG48           |SRR1027475    |PRJNA226980 |female |midgut         |48              |12_48                 |       61002|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_midgut         |NA                               |
|SG120          |SRR1027474    |PRJNA226980 |female |salivary_gland |120             |72_144                |       67277|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_salivary_gland |female_x_salivary_gland_x_72_144 |
|SG96           |SRR1027473    |PRJNA226980 |female |salivary_gland |96              |72_144                |       86616|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_salivary_gland |female_x_salivary_gland_x_72_144 |
|SG48           |SRR1027471    |PRJNA226980 |female |salivary_gland |48              |12_48                 |      179494|https://doi.org/10.1186/s13071-017-2080-1; https://doi.org/10.1186/1471-2164-15-518 |female_x_salivary_gland |female_x_salivary_gland_x_12_48  |
|NEm            |SRR4416251    |PRJNA327120 |male   |whole          |168             |168_264               |    39445131|https://doi.org/10.1016/j.aspen.2018.05.009                                         |male_x_whole            |NA                               |
|NEf            |SRR4416250    |PRJNA327120 |female |whole          |168             |168_264               |    28403332|https://doi.org/10.1016/j.aspen.2018.05.009                                         |female_x_whole          |NA                               |

#### Combining variables allows us to build two differential expression models

#### Exploring differential expression results using a Shiny App

### Lessons learned

#### Lessons learned from gene counting

### Next Steps

## About the GitHub repository 

### Repository Structure

- `inputs/`: Directory containing metadata and other necessary inputs for the workflow.
- `scripts/`: Directory containing scripts for data processing and analysis.
- `envs/`: Directory containing conda environment files for software dependencies.
- `Snakefile`: The main workflow script orchestrating the analysis.
- `notebooks`: Directory containing analysis notebooks. These notebooks are scratch analyses that were used to determine the best next steps for the Snakemake pipeline or for visualizations to include in the [Shiny App](./shiny).
- `shiny`: Directory containing the code and data required to run the *A. americanum* Differential Expression Explorer Shiny App. Also contains the instructions for running the app locally or via the internet. 

### Workflow

The main analysis is orchestrated by a `Snakefile` which defines various rules for processing the data.
Here is a brief description of the workflow:

1. **Download Data**: Downloads the raw sequencing data from SRA using `fasterq-dump`.
2. **Pre-processing**: Combines multiple run accessions into one file per Illumina library name, performs quality control and trimming with `fastp`, and separates paired-end reads.
3. **Read Quantification**: Indexes the transcriptome using `Salmon`, quantifies transcripts, and creates a tx2gene file by mapping transcripts back to the genome using `uLTRA` and assigning transcripts gene names by determining the gene each transcript overlaps.
4. **Differential Expression Analysis**: Builds differential expression models to study the effects of sex, tissue, and blood meal hour on gene expression. These differential expression models are then used as input to the [Shiny App](./shiny).

The [Snakefile](./Snakefile) contains additional documentation detailing the purpose of each rule.

### Running the pipeline on an AWS EC2 instance

We used the Canonical, Ubuntu, 22.04 LTS, amd64 jammy image build on 2023-05-16 with 64 bit architecture and AMI ID ami-0f8e81a3da6e2510a.
We initially launched an m5a.large instance, and after configuration ran the pipeline on m5a.2xlarge instance type.
To set up the instance to run the pipeline, we installed and configured miniconda with the commands below.

```
curl -JLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download the miniconda installation script
bash Miniconda3-latest-Linux-x86_64.sh # run the miniconda installation script. Accept the license and follow the defaults.
source ~/.bashrc # source the .bashrc for miniconda to be available in the environment

# configure miniconda channel order
conda config --add channels defaults 
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install mamba # install mamba for faster software installation.
```

### Running the snakemake pipeline

First, clone the repository to your computer.

```
git clone git@github.com:Arcadia-Science/2023-amblyomma-americanum-diffex.git # assumes SSH is configured
```

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html) (see above for linux).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.
We installaed Miniconda3 version `py311_23.5.2-0` and mamba version `1.4.9`.

```
mamba env create -n amam --file environment.yml
conda activate amam
```

To start the pipeline, run:

```
snakemake --use-conda -j 2
```

The pipeline processes accessions listed in a TSV file, which defaults to `inputs/metadata.tsv`.
This can be changed on line 4 of the [Snakefile](./Snakefile).

## Contributing

See [this guide](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md) to see how we recognize feedback and contributions on our code.

## IP decision and pub

Currently, we have been advised to keep this GitHub repo private and to not pub out this project.
This is because the differential expression analysis we performed is trivial and does not itself produce interesting and shareable biological insights, while the repo houses a [Shiny App](./shiny) that may be considered an asset for Trove.
The IP guidance may change over time.
In the mean time, the READMEs in this repository serve as the main documentation for the project, alongside the Notion Roadmaps in the Software Documentation section.

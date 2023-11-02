import pandas as pd

# read in metadata file
metadata_all = pd.read_csv("inputs/metadata.tsv", sep = "\t").set_index("run_accession", drop = False)
# select columns that we need metadata from for wildcards and other places in the workflow
metadata_illumina = metadata_all[["library_name", "condition", "library_layout"]].drop_duplicates()
# set the index to library name to allow dictionary-like lookups from the metadata tables with param lambda functions
metadata_illumina = metadata_illumina.set_index("library_name", drop = False)

# use metadata tables to create global variables
# extract SRA accessions to a variable, which we'll use to download the raw data
RUN_ACCESSIONS = metadata_all["run_accession"].unique().tolist()
# extract library names; some libraries are split between multiple SRA accessions
ILLUMINA_LIB_NAMES = metadata_illumina["library_name"].unique().tolist()

rule all:
    input: "outputs/tx2gene/tx2gene.tsv"

######################################
# Download short & long read data
######################################

rule download_fastq_files:
    output: temp("inputs/raw/{run_accession}.fq.gz")
    conda: "envs/sratools.yml"
    params: outdir = "inputs/raw/"
    shell:'''
    fasterq-dump --split-spot -Z {wildcards.run_accession} | gzip > {output}
    ''' 

######################################
## Process & assemble illumina files
######################################

rule combine_by_library_name:
    """
    Some of the input sequences have multiple run accessions (SRR*) associated with a single library.
    This rule combines those run accessions into one file by ill_lib_name (illumina library name).
    Since it uses the metadata_illumina file to do this, as well as expanding over the runacc wild card in the input,
    the output wildcard illumina_lib_name will only contain illumina library names.
    """
    input: expand("inputs/raw/{run_accession}.fq.gz", run_accession = RUN_ACCESSIONS)
    output: temp(expand("outputs/read_qc/raw_combined/{illumina_lib_name}.fq.gz", illumina_lib_name = ILLUMINA_LIB_NAMES))
    params: 
        indir = "inputs/raw/",
        outdir = "outputs/raw_combined/"
    run:
        # create a dictionary that has library names as keys and run accessions as values
        tmp_illumina = metadata_all[metadata_all["instrument"] != "Sequel II"].drop_duplicates()
        tmp = tmp_illumina[["library_name", "run_accession"]]
        libdict = {}
        for group, d in tmp.groupby('library_name'):
            libdict[group] = d['run_accession'].values.tolist()

        # format the values to be a string of file paths, each separated by a space
        for library_name, run_accessions in libdict.items():
            library_paths = []
            for run_accession in run_accessions:
                path = params.indir + run_accession + ".fq.gz"
                library_paths.append(path)
                shell_drop_in = " ".join(library_paths)
            shell("cat {shell_drop_in} > {params.outdir}/{library_name}.fq.gz")

rule fastp:
    """
    We set the quality trimming parameters used by the Oyster River Protocol for de novo transcriptomics:
    - quality trim with a phred score of 2
    - trim the polyA tails
    - adapter trim
    """
    input: "outputs/read_qc/raw_combined/{illumina_lib_name}.fq.gz"
    output:
        json = "outputs/read_qc/fastp/{illumina_lib_name}.json",
        html = "outputs/read_qc/fastp/{illumina_lib_name}.html",
        fq = "outputs/read_qc/fastp/{illumina_lib_name}.fq.gz"
    conda: "envs/fastp.yml"
    threads: 2
    params: liblayout = lambda wildcards: metadata_illumina.loc[wildcards.illumina_lib_name, "library_layout"]
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        fastp -i {input} --thread {threads} --trim_poly_x --qualified_quality_phred 2 --json {output.json} --html {output.html} --report_title {wildcards.illumina_lib_name} --interleaved_in --stdout | gzip > {output.fq}
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        fastp -i {input} --thread {threads} --trim_poly_x --qualified_quality_phred 2 --json {output.json} --html {output.html} --report_title {wildcards.illumina_lib_name} --stdout | gzip > {output.fq}
    fi
    '''

rule split_paired_end_reads_fastp:
    input: fq = "outputs/read_qc/fastp/{illumina_lib_name}.fq.gz"
    output:
        r1="outputs/read_qc/fastp_separated_reads/{illumina_lib_name}_R1.fq.gz",
        r2="outputs/read_qc/fastp_separated_reads/{illumina_lib_name}_R2.fq.gz"
    conda: "envs/bbmap.yml"
    params: liblayout = lambda wildcards: metadata_illumina.loc[wildcards.illumina_lib_name, "library_layout"]
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        repair.sh in={input} out={output.r1} out2={output.r2} repair=t overwrite=true -Xmx60g
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        cp {input} {output.r1}
        touch {output.r2}
    fi
    '''

##############################################
## Read quantification
##############################################

# using EVM genome & genome annotation.
# for now, downloaded by hand from S3 until available on zenodo.
# documented for myself in a README on AWS in the inputs folder.

# rule download_transcriptome:
# will add when there's a public download link

rule index_transcriptome:
    input: "inputs/assembly/orthofuser_final_clean.fa.dammit.fasta"
    output: "outputs/quantification/salmon_index/info.json"
    threads: 1
    params: indexdir = "outputs/quantification/salmon_index/"
    conda: "envs/salmon.yml"
    shell:'''
    salmon index -t {input} -i {params.indexdir} -k 31
    '''

rule salmon:
    input:
        index = "outputs/quantification/salmon_index/info.json",
        r1="outputs/read_qc/fastp_separated_reads/{illumina_lib_name}_R1.fq.gz",
        r2="outputs/read_qc/fastp_separated_reads/{illumina_lib_name}_R2.fq.gz"
    output: 
        "outputs/quantification/salmon/{illumina_lib_name}_quant/quant.sf",
        "outputs/quantification/salmon/{illumina_lib_name}_quant/aux_info/eq_classes.txt.gz"
    params: 
        liblayout = lambda wildcards: metadata_illumina.loc[wildcards.illumina_lib_name, "library_layout"],
        indexdir = "outputs/quantification/salmon_index/",
        outdir = lambda wildcards: "outputs/quantification/salmon/" + wildcards.illumina_lib_name + "_quant" 
    conda: "envs/salmon.yml"
    threads: 2
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        salmon quant -i {params.indexdir} -l A -1 {input.r1} -2 {input.r2} -o {params.outdir} --dumpEq --writeOrphanLinks -p {threads} 
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        salmon quant -i {params.indexdir} -l A -r {input.r1} -o {params.outdir} --dumpEq --writeOrphanLinks -p {threads}
    fi
    '''

# create a tx2gene file by mapping transcripts back to genome with a splice-aware long read aligner

rule convert_gff_to_gtf:
    input: "inputs/annotations/evm/Amblyomma_americanum_filtered_assembly.evm.gff3"
    output: "inputs/annotations/evm/Amblyomma_americanum_filtered_assembly.evm.gtf"
    conda: "envs/agat.yml"
    shell:'''
    agat_convert_sp_gff2gtf.pl --gff {input} -o {output}
    '''

rule map_transcripts_to_genome_with_ultra:
    input:
        genome="inputs/genome/Amblyomma_americanum_filtered_assembly.fasta",
        gtf="inputs/annotations/evm/Amblyomma_americanum_filtered_assembly.evm.gtf",
        txome="inputs/assembly/orthofuser_final_clean.fa.dammit.fasta"
    output: "outputs/tx2gene/ultra/Amblyomma_americanum_filtered_assembly.sam"
    params: 
        outdir = "outputs/tx2gene/ultra/",
        prefix = "Amblyomma_americanum_filtered_assembly"
    threads: 14
    conda: "envs/ultra.yml"
    shell:'''
    uLTRA pipeline {input.genome} {input.gtf} {input.txome} {params.outdir} --prefix {params.prefix} --t 14 --disable_infer
    '''

rule assign_transcripts_to_genes_by_overlaps_with_gtf_genes:
    input: 
        sam="outputs/tx2gene/ultra/Amblyomma_americanum_filtered_assembly.sam",
        gtf="inputs/annotations/evm/Amblyomma_americanum_filtered_assembly.evm.gtf"
    output: "outputs/tx2gene/tx2gene.tsv"
    conda: "envs/pysam.yml"
    shell:'''
    python scripts/assign_mapped_transcripts_to_gene_by_gtf_overlap.py {input.gtf} {input.sam} {output}
    '''

rule make_counts:

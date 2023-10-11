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
    input: "outputs/counts/raw_counts.tsv"

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

# using EVM genome annotation.
# for now, downloaded by hand from S3 until available on zenodo.

rule star_index_genome:
    input:
        genome = "inputs/genome/Amblyomma_americanum_filtered_assembly.fasta",
        gff = "inputs/annotations/evm/Amblyomma_americanum_filtered_assembly.evm.gff3" 
    output: "inputs/genome/SAindex"
    conda: "envs/star.yml"
    threads: 16
    shell:'''
    STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir genome  \
         --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gff} \
         --sjdbGTFtagExonParentTranscript mRNA --sjdbOverhang  99
    '''

rule star_map_reads:
    input:
        genome_index = "inputs/genome/SAindex",
        fastq = expand("outputs/read_qc/fastp_separated_reads/{{illumina_lib_name}}_R{pair}.fq.gz", pair = [1, 2])
    output: "outputs/counts/star/{illumina_lib_name}_Aligned.sortedByCoord.out.bam"
    params: 
        genomedir = "inputs/genome",
        outprefix = lambda wildcards: "outputs/counts/star/" + wildcards.illumina_lib_name + "_",
        liblayout = lambda wildcards: metadata_illumina.loc[wildcards.illumina_lib_name, "library_layout"]
    threads: 4
    conda: "envs/star.yml"
    shell:'''
    # define a bash variable so the STAR command only needs to be repeated once
    if [ "{params.liblayout}" == "PAIRED" ]; then
        fastq_files={input.fastq}
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        fastq_files={input.fastq[0]}
    fi

    STAR --runThreadN {threads} {params.genomedir}            \
         --readFilesIn ${{fastq_files}} --outFilterType BySJout  \
         --outFilterMultimapNmax 20 --alignSJoverhangMin 8    \
         --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 \
         --alignIntronMax 1000000 --alignMatesGapMax 1000000  \
         --outSAMattributes NH HI NM MD --outSAMtype BAM      \
         SortedByCoordinate --outFileNamePrefix {params.outprefix}
    '''

rule samtools_index:
    input:"outputs/counts/star/{illumina_lib_name}_Aligned.sortedByCoord.out.bam"
    output:"outputs/counts/star/{illumina_lib_name}_Aligned.sortedByCoord.out.bam.bai"
    conda: "envs/samtools.yml"
    shell:'''
    samtools index {input}
    '''

rule htseq_count:
    input:
        bai="outputs/counts/star/{illumina_lib_name}_Aligned.sortedByCoord.out.bam.bai",
        gff="inputs/annotations/evm/Amblyomma_americanum_filtered_assembly.evm.gff3"
    output: "outputs/counts/htseq_count/{illumina_lib_name}.out"
    conda: "envs/htseq.yml"
    shell:'''
    htseq-count -f bam {input.bai} {input.gff} -i Parent -r pos > {output}
    '''

rule combine_htseq_counts:
    input: expand("outputs/counts/htseq_count/{illumina_lib_name}.out", illumina_lib_name = ILLUMINA_LIB_NAMES)
    output: "outputs/counts/raw_counts.tsv"
    conda: "envs/tidyverse.yml"
    shell:'''
    Rscript bin/combine_htseq_counts.R {output} {input}
    '''

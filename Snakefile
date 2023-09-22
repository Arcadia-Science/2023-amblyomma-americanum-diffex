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
    input: "outputs/quantification/grouper/mag.flat.clust"

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

rule download_transcriptome:

rule index_transcriptome:
    input: "inputs/assembly/orthofuser_final_clean.fa.dammit.fasta"
    output: "outputs/quantification/salmon_index/info.json"
    threads: 1
    params: indexdir = "outputs/quantification/salmon_index/"
    conda: "envs/salmon.yml"
    shell:'''
    salmon index -t {input} -i {params.indexdir} -k 31
    '''

rule salmon_for_grouper:
    input:
        index = "outputs/quantification/salmon_index/info.json",
        r1="outputs/read_qc/fastp_separated_reads/{illumina_lib_name}_R1.fq.gz",
        r2="outputs/read_qc/fastp_separated_reads/{illumina_lib_name}_R2.fq.gz"
    output: "outputs/quantification/salmon/{illumina_lib_name}_quant/quant.sf"
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

rule make_grouper_config_file:
    """
    Grouper requires a config file in the following format:
    conditions:
        - Control
        - HOXA1 Knockdown
    samples:
        Control:
            - SRR493366_quant
            - SRR493367_quant
        HOXA1 Knockdown:
            - SRR493369_quant
            - SRR493371_quant
    outdir: human_grouper
    orphan: True
    mincut: True
    """
    input: expand("outputs/quantification/salmon/{illumina_lib_name}_quant/quant.sf", illumina_lib_name = ILLUMINA_LIB_NAMES)
    output: conf = "outputs/quantification/grouper/grouper_conf.yml"
    params: 
        grouperdir = "outputs/quantification/grouper/",
        salmondir =  "outputs/quantification/salmon/"
    run:
        # create a dictionary of assembly groups: library names
        tmp = metadata_illumina[["condition", "library_name"]]
        condition_dict = {}
        for group, d in tmp.groupby('condition'):
            condition_dict[group] = d['library_name'].values.tolist()
        # use the dictionary to parse a string of conditions (assembly groups) that will be written to the grouper yaml
        conditions_list = "\n    - ".join(list(condition_dict.keys()))
        # use the dictionary to parse a nested string of conditions: salmon results by library name that will be written to grouper yaml
        samples_list = []
        for condition, library_names in condition_dict.items():
            samples_list.append("\n    - " + condition)
            for library_name in library_names:
                samples_list.append("\n        -" + params.salmondir + library_name + "_quant")

        samples_list = "".join(samples_list)
        # create a config template with format strings that will be substituted in the write process
        config_template = """\
conditions:
    - {conditions_list}
samples: {samples_list}
outdir: {outdir}
orphan: True
mincut: True
"""
        with open(output.conf, 'wt') as fp:
            fp.write(config_template.format(conditions_list = conditions_list, 
                                            samples_list = samples_list,
                                            outdir = params.grouperdir))


rule run_grouper:
    input: "outputs/quantification/grouper/grouper_conf.yml"
    output: "outputs/quantification/grouper/mag.flat.clust"
    conda: "envs/biogrouper.yml"
    shell:'''
    Grouper --config {input}
    '''

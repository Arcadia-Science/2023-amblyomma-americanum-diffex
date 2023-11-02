#!/usr/bin/env python3c

import pysam
import pyranges as pr
import argparse

def gtf_to_pyranges(gtf_file):
    """Convert a GTF file to a PyRanges object."""
    df = pr.read_gtf(gtf_file).df  # Retrieve the DataFrame
    genes = df[df["Feature"] == "gene"].copy()
    genes["gene_name"] = genes["gene_id"].astype(str)
    return pr.PyRanges(genes[["Chromosome", "Start", "End", "Strand", "gene_name"]])

def sam_to_pyranges(sam_file):
    """Convert a SAM/BAM file to a PyRanges object."""
    with pysam.AlignmentFile(sam_file, "r") as sam:
        data = []
        for read in sam:
            if not read.is_unmapped:
                data.append([read.reference_name, read.reference_start, read.reference_end, read.query_name])
        df = pr.pd.DataFrame(data, columns=["Chromosome", "Start", "End", "Name"])
        return pr.PyRanges(df)


def calculate_overlap(Start, End, Start_gene, End_gene, gene_name):
    no_overlap = gene_name == "-1"
    overlap_start = pr.pd.Series([max(s, g) for s, g in zip(Start, Start_gene)])
    overlap_end = pr.pd.Series([min(e, g) for e, g in zip(End, End_gene)])
    
    overlap = overlap_end - overlap_start
    overlap[overlap < 0] = 0
    overlap[no_overlap] = 0

    return overlap

# test case:
# test_row = {'gene_name': 'evm.TU.contig_177737_1.1', 'Start': 0, 'End': 14685, 'Start_gene': 5972, 'End_gene': 27724}
# print(calculate_overlap(test_row))

def main(gtf_file, sam_file, output_file):
    gene_ranges = gtf_to_pyranges(gtf_file)
    transcript_ranges = sam_to_pyranges(sam_file)

    # Joining with 'outer' will keep transcripts without gene overlap as well
    overlaps = transcript_ranges.join(gene_ranges, how='outer', suffix="_gene")
    overlaps_df = overlaps.df

    # Calculate overlap length
    overlaps_df['overlap_length'] = calculate_overlap(
        overlaps_df['Start'], overlaps_df['End'], 
        overlaps_df['Start_gene'], overlaps_df['End_gene'], 
        overlaps_df['gene_name']
    )

    print("Columns in overlaps dataframe:", overlaps_df.columns)

    # Keep only the gene with the largest overlap for each transcript
    overlaps_df = overlaps_df.loc[overlaps_df.groupby('Name')['overlap_length'].idxmax()]


    with open(output_file, 'w') as out:
        for index, row in overlaps_df.iterrows():
            # If there's no gene overlap, use the transcript name as gene name
            #gene_name = row['gene_name'] if row['gene_name'] != "-1" else row['Name']
            #out.write(f"{row['Name']}\t{gene_name}\n")
            # write only transcripts that overlapped with a gene
            out.write(f"{row['Name']}\t{row['gene_name']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map transcripts to genes based on overlap.")
    parser.add_argument("gtf_file", help="Path to the GTF annotation file.")
    parser.add_argument("sam_file", help="Path to the SAM/BAM file.")
    parser.add_argument("output_file", help="Path to the output file.")
    args = parser.parse_args()

    main(args.gtf_file, args.sam_file, args.output_file)

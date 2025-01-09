#!/usr/bin/env python3

import os
import pandas as pd
from Bio import SeqIO

def get_genome_length(fasta_file):
    """Calculate the total genome length from a reference FASTA file."""
    genome_length = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))
    print(f"Total genome length: {genome_length}")
    return genome_length

def calculate_metrics(input_file, genome_length):
    """Calculates the requested metrics based on the input TSV file."""
    # Load the TSV file into a DataFrame
    df = pd.read_csv(input_file, sep='\t')
    
    # Compute metrics
    results = []
    for _, row in df.iterrows():
        # Extract relevant fields
        sample_id = row.get("sample_id", "")
        total_length = row.get("total_length_pf", 0)
        bases_mapped = row.get("bases_mapped_pf", 0)
        mismatches = row.get("mismatches_pf", 0)
        
        # Calculations based on WDL logic
        aligned_frac_bases = bases_mapped / total_length if total_length != 0 else 0.0
        aligned_est_fold_cov = bases_mapped / genome_length if genome_length != 0 else 0.0
        average_identity = (
            100.0 - (100.0 * mismatches / bases_mapped) if bases_mapped != 0 else 0.0
        )
        
        # Add calculated values to results
        results.append({
            "SampleID": sample_id,
            "aligned_frac_bases": aligned_frac_bases,
            "aligned_est_fold_cov": aligned_est_fold_cov,
            "average_identity": average_identity,
            "total_length": total_length,
            "bases_mapped": bases_mapped,
            "mismatches": mismatches
        })
    
    return results

def main(base_dir, fasta_file):
    """Search directories ending with QC_results and process their files."""
    # Get genome length from FASTA file
    genome_length = get_genome_length(fasta_file)
    
    output_results = []

    for root, dirs, files in os.walk(base_dir):
        if root.endswith("QC_results"):
            file_path = os.path.join(root, "final_qc_reports/Bam_stats_pf_Final.tsv")
            if os.path.exists(file_path):
                print(f"Processing {file_path}...")
                metrics = calculate_metrics(file_path, genome_length)
                output_results.extend(metrics)
    
    # Create a DataFrame for results
    results_df = pd.DataFrame(output_results)
    
    # Save the results to a CSV file
    output_file = "KaramokoQC.csv"
    results_df.to_csv(output_file, index=False)
    print(f"Metrics summary saved to {output_file}")


if __name__ == "__main__":
    # Define the base directory and reference genome path
    base_directory = "/data/general/finterly/S3TIMULATE_WGS"
    reference_fasta = "PlasmoDB-61_Pfalciparum3D7_Genome.fasta"
    main(base_directory, reference_fasta)

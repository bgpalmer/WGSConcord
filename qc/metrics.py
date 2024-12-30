#!/usr/bin/env python3

import os
import pandas as pd

def calculate_metrics(input_file):
    """Calculates the requested metrics based on the input TSV file."""
    # Load the TSV file into a DataFrame
    df = pd.read_csv(input_file, sep='\t')
    
    # Compute metrics
    results = []
    for _, row in df.iterrows():
        genome_length = row.get("total_length_pf", 0)
        bases_mapped = row.get("bases_mapped_pf", 0)
        total_bases = row.get("bases_mapped_(cigar)_pf", 0)
        mismatches = row.get("mismatches_pf", 0)
        
        # Raw estimated fold coverage
        raw_est_fold_cov = total_bases / genome_length if genome_length != 0 else 0.0
        
        # Aligned fraction of bases
        aligned_frac_bases = bases_mapped / total_bases if total_bases != 0 else 0.0
        
        # Aligned estimated fold coverage
        aligned_est_fold_cov = bases_mapped / genome_length if genome_length != 0 else 0.0
        
        # Average identity
        average_identity = (
            100.0 - (100.0 * mismatches / bases_mapped)
            if bases_mapped != 0
            else 0.0
        )
        
        results.append({
            "sample_id": row["sample_id"],
            "raw_est_fold_cov": raw_est_fold_cov,
            "aligned_frac_bases": aligned_frac_bases,
            "aligned_est_fold_cov": aligned_est_fold_cov,
            "average_identity": average_identity,
        })
    
    return results


def main(base_dir):
    """Search directories ending with QC_results and process their files."""
    output_results = []

    for root, dirs, files in os.walk(base_dir):
        if root.endswith("QC_results"):
            file_path = os.path.join(root, "Bam_stats_pf_Final.tsv")
            if os.path.exists(file_path):
                print(f"Processing {file_path}...")
                metrics = calculate_metrics(file_path)
                output_results.extend(metrics)
    
    # Create a DataFrame for results
    results_df = pd.DataFrame(output_results)
    
    # Save the results to a CSV file
    output_file = os.path.join(base_dir, "QC_metrics_summary.csv")
    results_df.to_csv(output_file, index=False)
    print(f"Metrics summary saved to {output_file}")


if __name__ == "__main__":
    # Define the base directory (adjust to your path)
    base_directory = "/data/general/finterly/S3TIMULATE_WGS"
    main(base_directory)

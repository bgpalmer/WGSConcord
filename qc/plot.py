#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import re

def load_data(file_path):
    """Load data from a CSV file."""
    return pd.read_csv(file_path)

def find_common_samples(df1, df2):
    """Find common samples between two DataFrames based on 'SampleID'."""
    return pd.merge(df1, df2, on="SampleID", suffixes=("_Terra", "_Karamoko"))

def extract_batch_and_sample(sample_id):
    """Extract batch and sample number from SampleID."""
    batch_match = re.search(r"batch[A-Z]", sample_id)
    sample_match = re.search(r"S\d+$", sample_id)
    
    batch = batch_match.group(0) if batch_match else "Unknown"
    sample = sample_match.group(0) if sample_match else "Unknown"
    return batch, sample

def plot_metrics(common_df, metrics):
    """Plot metrics for common samples."""
    for metric in metrics:
        terra_metric = f"{metric}_Terra"
        karamoko_metric = f"{metric}_Karamoko"
        
        if terra_metric in common_df.columns and karamoko_metric in common_df.columns:
            # Determine axis limits to ensure same scale
            min_value = min(common_df[terra_metric].min(), common_df[karamoko_metric].min())
            max_value = max(common_df[terra_metric].max(), common_df[karamoko_metric].max())
            
            plt.figure(figsize=(10, 8))
            plt.scatter(common_df[terra_metric], common_df[karamoko_metric], alpha=0.7)
            plt.plot([min_value, max_value], [min_value, max_value], 'r--', label="Perfect Agreement")
            
            # Annotate each point with batch and sample
            for _, row in common_df.iterrows():
                batch, sample = extract_batch_and_sample(row["SampleID"])
                plt.text(
                    row[terra_metric],
                    row[karamoko_metric],
                    f"{batch}-{sample}",
                    fontsize=8,
                    alpha=0.7
                )
            
            # Set the axis limits
            plt.xlim(min_value, max_value)
            plt.ylim(min_value, max_value)
            
            plt.title(f"{metric} Comparison")
            plt.xlabel(f"Terra {metric}")
            plt.ylabel(f"Karamoko {metric}")
            plt.legend()
            plt.grid(True)
            plt.savefig(f"{metric}_comparison.png")
            plt.show()

def main(terra_file, karamoko_file):
    """Compare metrics between Terra and Karamoko."""
    # Load data
    terra_df = load_data(terra_file)
    karamoko_df = load_data(karamoko_file)
    
    # Find common samples
    common_df = find_common_samples(terra_df, karamoko_df)
    print(f"Found {len(common_df)} common samples.")
    
    # Define metrics to compare
    metrics = [
        "raw_est_fold_cov",
        "aligned_frac_bases",
        "aligned_est_fold_cov",
        "average_identity"
    ]
    
    # Plot metrics
    plot_metrics(common_df, metrics)

if __name__ == "__main__":
    # Update these paths to the actual CSV file locations
    terra_csv = "TerraQC.csv"
    karamoko_csv = "KaramokoQC.csv"
    main(terra_csv, karamoko_csv)

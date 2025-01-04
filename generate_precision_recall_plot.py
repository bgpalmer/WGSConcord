import sys
import os
import re
import matplotlib.pyplot as plt

def parse_stats_file(stats_file):
    """Parse a bcftools stats file for TP, FP, and FN values."""
    tp, fp, fn = 0, 0, 0
    records = {}
    with open(stats_file, 'r') as f:
        for line in f:
            if line.startswith("SN"):
                parts = line.strip().split("\t")
                if len(parts) >= 4 and "number of records" in parts[2]:
                    set_id = int(parts[1])
                    records[set_id] = int(parts[3])

    print(records)
    tp = records.get(2, 0)
    fp = records.get(0, 0) - tp
    fn = records.get(1, 0) - tp
    return tp, fp, fn

def calculate_precision_recall(tp, fp, fn):
    """Calculate precision and recall."""
    precision = tp / (tp + fp) if tp + fp > 0 else 0
    recall = tp / (tp + fn) if tp + fn > 0 else 0
    return precision, recall

def generate_plot(sample_key, chromosomes, precisions, recalls, output_plot):
    """Generate precision and recall plot."""
    plt.figure(figsize=(10, 6))
    plt.plot(chromosomes, precisions, label="Precision", marker="o", linestyle="--")
    plt.plot(chromosomes, recalls, label="Recall", marker="o", linestyle="-")
    plt.xlabel("Chromosome")
    plt.ylabel("Metric Value")
    plt.title(f"Precision and Recall by Chromosome for {sample_key}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()

def main(sample_key, stats_dir, output_plot):
    chromosomes = []
    precisions = []
    recalls = []

    for stats_file in sorted(os.listdir(stats_dir)):
        if stats_file.startswith(sample_key) and stats_file.endswith(".stats.txt"):
            chrom_match = re.search(r"_Pf3D7_(\d+)_v3", stats_file)
            if chrom_match:
                chrom = chrom_match.group(1)
                tp, fp, fn = parse_stats_file(os.path.join(stats_dir, stats_file))
                print(tp, fp, fn)
                precision, recall = calculate_precision_recall(tp, fp, fn)
                chromosomes.append(f"Chr{chrom}")
                precisions.append(precision)
                recalls.append(recall)

    if chromosomes:
        generate_plot(sample_key, chromosomes, precisions, recalls, output_plot)
        print(f"Precision and recall plot saved for sample {sample_key}.")
    else:
        print(f"No valid statistics found for sample {sample_key}.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 generate_precision_recall_plot.py <sample_key> <stats_dir> <output_plot>")
        sys.exit(1)

    sample_key = sys.argv[1]
    stats_dir = sys.argv[2]
    output_plot = sys.argv[3]
    main(sample_key, stats_dir, output_plot)

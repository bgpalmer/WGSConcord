#!/usr/bin/env python3

import cairo
import os
import pandas as pd
import re
import logging
from math import ceil
import colorsys
import matplotlib.pyplot as plt

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("visualization_debug.log"),
        logging.StreamHandler()
    ]
)


def get_chromosome_lengths(ref_genome_path):
    """
    Parse the reference genome to get chromosome lengths.
    """
    lengths = {}
    try:
        with open(ref_genome_path, "r") as f:
            chrom = None
            bases = 0
            for line in f:
                if line.startswith(">"):
                    if chrom:
                        lengths[chrom] = bases
                        logging.info(f"Chromosome {chrom} length: {bases} bases")
                    chrom = re.search(r"Pf3D7_(\d+)", line)
                    chrom = chrom.group(1) if chrom else None
                    bases = 0
                else:
                    bases += len(line.strip())
            if chrom:
                lengths[chrom] = bases
                logging.info(f"Chromosome {chrom} length: {bases} bases")
    except Exception as e:
        logging.error(f"Error reading reference genome: {e}")
    return lengths


def parse_concordance_files(concordance_dir):
    """
    Parse all snpsift_concordance.txt files and organize by chromosome.
    """
    data = {}
    samples = set()
    try:
        for root, _, files in os.walk(concordance_dir):
            for file in files:
                if file.endswith("snpsift_concordance.txt"):
                    sample = os.path.basename(root)
                    
                    # Extract chromosome
                    chromosome_match = re.search(r"_Pf3D7_(\d+)", sample)
                    if not chromosome_match:
                        logging.warning(f"Could not determine chromosome from sample: {sample}")
                        continue
                    chromosome = chromosome_match.group(1)

                    # Extract Sample ID (S#)
                    sample_id_match = re.search(r"(S\d+)_Pf3D7_", sample)
                    batch_id_match = re.search(r"(batch[A-Z])", sample)
                    if not sample_id_match:
                        logging.warning(f"Could not determine sample ID from: {sample}")
                        continue
                    sample_id = batch_id_match.group(1) + "_" + sample_id_match.group(1)
                    print(sample_id)
                    samples.add(sample_id)

                    file_path = os.path.join(root, file)
                    logging.info(f"Processing file: {file_path}")

                    try:
                        # Read the file and handle header mismatch
                        with open(file_path, "r") as f:
                            header = f.readline().strip().split("\t")
                            df = pd.read_csv(f, sep="\t", names=header, dtype=str, index_col=False)

                        # Skip rows with '#Total' in the 'chr' column
                        if "chr" in df.columns:
                            df = df[df["chr"] != "#Total"]

                        if "pos" not in df.columns:
                            logging.warning(f"'pos' column missing in file: {file_path}")
                            continue

                        # Convert the 'pos' column to numeric and handle invalid positions
                        df["Position"] = pd.to_numeric(df["pos"], errors="coerce")
                        invalid_positions = df[pd.to_numeric(df["pos"], errors="coerce").isna()]
                        if not invalid_positions.empty:
                            logging.warning(
                                f"File {file_path} contains {len(invalid_positions)} invalid positions. "
                                f"Examples: {invalid_positions[['pos']].head().to_dict(orient='records')}"
                            )
                        df.dropna(subset=["Position"], inplace=True)

                        # Add extracted chromosome and sample ID
                        df["Chromosome"] = chromosome
                        df["Sample"] = sample_id

                        if chromosome not in data:
                            data[chromosome] = []
                        data[chromosome].append(df)
                    except Exception as e:
                        logging.error(f"Error processing file {file_path}: {e}")
    except Exception as e:
        logging.error(f"Error parsing concordance files: {e}")

    # Limit to 10 unique samples
    samples = sorted(samples)
    logging.info(f"Selected samples: {samples}")
    return data, samples


def create_combined_chromosome_page(chromosome_data, chromosome_lengths, samples, output_file):
    """
    Create a single-page visualization with a dark theme and vibrant color-blind-friendly palette.
    """
    logging.info("Creating combined visualization with dark theme and accessible colors")
    try:
        # Canvas dimensions
        num_chromosomes = len(chromosome_lengths)
        num_samples = len(samples)
        width = 2000
        height = max(200 * num_chromosomes + 300, 100 * num_samples + 300)
        margin = 150
        max_length = max(chromosome_lengths.values())
        x_scale = (width - 2 * margin) / max_length
        y_spacing = 200  # Space between chromosomes
        sample_offset = 10  # Offset for sample annotations

        logging.info(f"Canvas dimensions: width={width}, height={height}, margin={margin}")
        logging.info(f"Max chromosome length: {max_length}, X scale: {x_scale:.2f}")

        # Create a surface and context
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        ctx = cairo.Context(surface)

        # Set dark background
        ctx.rectangle(0, 0, width, height)
        ctx.set_source_rgb(0.1, 0.1, 0.1)  # Dark gray background
        ctx.fill()
        logging.info("Canvas background set to dark gray")

        # Title
        ctx.set_source_rgb(1, 1, 1)  # White text
        ctx.set_font_size(28)
        ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        ctx.move_to(margin, 50)
        ctx.show_text("Chromosome Visualization with Sample Annotations")
        logging.info("Title added to canvas")

        # Use Matplotlib's tableau-colorblind palette
        cmap = plt.cm.get_cmap("tab10")  # Alternatively, use "tab20" for more colors
        sample_colors = {
            sample: cmap(i / max(1, len(samples)))[:3]  # Get RGB and convert to 0-1 scale
            for i, sample in enumerate(samples)
        }

        # Legend
        legend_x = width - 300
        legend_y = 100  # Initialize legend position explicitly
        ctx.set_source_rgb(1, 1, 1)  # White text for legend
        ctx.set_font_size(24)
        ctx.move_to(legend_x, legend_y)
        ctx.show_text("Legend:")

        for i, (sample_name, color) in enumerate(sample_colors.items()):
            ctx.set_source_rgb(*color)
            ctx.rectangle(legend_x + 20, legend_y + (i + 1) * 30 - 10, 20, 20)
            ctx.fill()
            ctx.set_source_rgb(1, 1, 1)  # White text for sample names
            ctx.move_to(legend_x + 50, legend_y + (i + 1) * 30 + 5)
            ctx.show_text(sample_name)
        logging.info("Legend added to canvas")

        # Plot each chromosome
        y_position = 150  # Initialize y_position for chromosomes
        ctx.move_to(margin - 50, 150)
        for idx, (chromosome, length) in enumerate(sorted(chromosome_lengths.items(), key=lambda x: int(x[0]))):
            logging.info(f"Processing Chromosome {chromosome}: Length={length}")

            # Chromosome label
            label_x = margin - 50  # Fixed x-coordinate for chromosome labels
            label_y = y_position  # Use current y_position for the label
            ctx.save()
            ctx.set_source_rgb(1, 1, 1)  # White text for chromosome labels
            ctx.set_font_size(24)
            ctx.translate(label_x, label_y)  # Set translation for the label
            ctx.rotate(-3.14159 / 2)  # Rotate 90 degrees for vertical label
            ctx.show_text(f"{chromosome}")
            ctx.restore()
            logging.info(f"Chromosome {chromosome} label added")

            # Annotations for each sample
            if chromosome in chromosome_data:
                samples_data = chromosome_data[chromosome]
                for sample_index, sample_name in enumerate(samples):  # Loop over predefined sample order
                    sample_data_list = [df for df in samples_data if df["Sample"].iloc[0] == sample_name]
                    if not sample_data_list:
                        continue
                    sample_data = sample_data_list[0]
                    color = sample_colors[sample_name]
                    ctx.set_source_rgb(*color)
                    ctx.set_line_width(0.15)

                    # Offset annotations for visibility
                    sample_y_offset = y_position + (sample_index - num_samples / 2) * sample_offset

                    for _, row in sample_data.iterrows():
                        pos = int(row["Position"])
                        x_pos = margin + pos * x_scale

                        # Draw smaller point
                        ctx.arc(x_pos, sample_y_offset, 0.15, 0, 2 * 3.14159)
                        ctx.fill()

            # Move to the next chromosome position
            y_position += y_spacing
            logging.info(f"Moving to next chromosome position: Y={y_position}")

        # Save to file
        surface.write_to_png(output_file)
        logging.info(f"Visualization saved to file: {output_file}")
    except Exception as e:
        logging.error(f"Error creating combined visualization: {e}")



# Main visualization function
def create_combined_visualization(chromosome_data, chromosome_lengths, samples, output_dir):
    """
    Create a combined visualization for all chromosomes on one page.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "combined_chromosome_visualization.png")
    create_combined_chromosome_page(chromosome_data, chromosome_lengths, samples, output_file)


def filter_concordance_data(concordance_data, mode="genotype_differences"):
    """
    Filter concordance data based on the mode.
    :param concordance_data: Parsed concordance data organized by chromosome.
    :param mode: Filter mode - "genotype_differences" or "missing_entries".
    :return: Filtered concordance data.
    """
    filtered_data = {}
    for chromosome, dfs in concordance_data.items():
        filtered_dfs = []
        for df in dfs:
            if mode == "genotype_differences":
                # Keep rows where both entries exist and genotypes differ
                genotype_columns = [col for col in df.columns if ("ALT" in col or "REF" in col) and "MISSING" not in col]
                filtered_df = df[df[genotype_columns].astype(int).sum(axis=1) > 0]  # Any genotype difference
            elif mode == "missing_entries":
                # Keep rows where one dataset is missing an entry
                missing_columns = [col for col in df.columns if "MISSING" in col]
                filtered_df = df[df[missing_columns].astype(int).sum(axis=1) > 0]  # Any missing entry
            else:
                logging.warning(f"Invalid mode: {mode}")
                filtered_df = pd.DataFrame()  # Empty DataFrame for invalid mode

            if not filtered_df.empty:
                filtered_dfs.append(filtered_df)

        if filtered_dfs:
            filtered_data[chromosome] = filtered_dfs

    return filtered_data


def create_combined_visualizations_by_mode(concordance_data, chromosome_lengths, samples, output_dir):
    """
    Create separate visualizations for genotype differences and missing entries.
    :param concordance_data: Parsed concordance data.
    :param chromosome_lengths: Chromosome lengths for scaling.
    :param samples: List of samples.
    :param output_dir: Directory to save the visualizations.
    """
    # Filter and visualize for genotype differences
    genotype_diff_data = filter_concordance_data(concordance_data, mode="genotype_differences")
    genotype_diff_file = os.path.join(output_dir, "genotype_differences_visualization.png")
    create_combined_chromosome_page(genotype_diff_data, chromosome_lengths, samples, genotype_diff_file)
    logging.info(f"Genotype differences visualization saved: {genotype_diff_file}")

    # Filter and visualize for missing entries
    missing_entries_data = filter_concordance_data(concordance_data, mode="missing_entries")
    missing_entries_file = os.path.join(output_dir, "missing_entries_visualization.png")
    create_combined_chromosome_page(missing_entries_data, chromosome_lengths, samples, missing_entries_file)
    logging.info(f"Missing entries visualization saved: {missing_entries_file}")


# Main
if __name__ == "__main__":
    # Configuration
    concordance_dir = "gvcf_comparison_results/concordance"
    reference_genome_path = "qc/PlasmoDB-61_Pfalciparum3D7_Genome.fasta"
    output_dir = "chromosome_visualizations"

    logging.info("Starting visualization script")

    # Parse reference genome and concordance data
    chromosome_lengths = get_chromosome_lengths(reference_genome_path)
    concordance_data, samples = parse_concordance_files(concordance_dir)

    # Create separate visualizations
    create_combined_visualizations_by_mode(concordance_data, chromosome_lengths, samples, output_dir)

    logging.info("Script completed successfully")

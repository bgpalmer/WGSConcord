#!/usr/bin/env python3

import cairo
import os
import pandas as pd
import re
import logging
import pysam
from math import ceil
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

def parse_isec_vcfs(isec_dir):
    """
    Parse VCF files from bcftools isec output directories.
    """
    data = {}
    samples = set()

    for root, _, files in os.walk(isec_dir):
        for file in files:
            if file.endswith("0000.vcf") or file.endswith("0001.vcf"):
                vcf_path = os.path.join(root, file)
                logging.info(f"Processing VCF: {vcf_path}")

                directory = os.path.dirname(vcf_path).split("/")[-1]
                print(directory)
                match = re.search(r"S3_WGS-(batch[A-Z]).*_(S\d+)_Pf3D7_(\d+)", directory)
                match2 = re.search(r"(S3_WGS-batch[A-Z].*_S\d+)_Pf3D7_.*", directory)
                if not match:
                    logging.warning(f"Could not extract sample or chromosome from {directory}")
                    continue

                batch, sn, chromosome = match.groups()
                sample = match2.groups()[0]

                # print(batch)
                # print(sn)
                # print(chromosome)
                # print(sample)

                sample_id = batch + "_" + sn
                samples.add(sample_id)
                if chromosome not in data:
                    data[chromosome] = []

                try:
                    vcf = pysam.VariantFile(vcf_path)
                    vcf_data = []
                    for record in vcf:
                        pos = record.pos
                        ref = record.ref
                        alts = record.alts if record.alts else []
                        gt = record.samples[sample]["GT"]
                        info = {
                            "Position": pos,
                            "REF": ref,
                            "ALT": ",".join(alts),
                            "Genotype": "/".join(map(str, gt)) if gt else ".",
                            "Sample": sample_id,
                            "Chromosome": chromosome
                        }
                        vcf_data.append(info)

                    df = pd.DataFrame(vcf_data)
                    data[chromosome].append(df)
                except Exception as e:
                    logging.error(f"Error parsing VCF {vcf_path}: {e}")

    return data, list(samples)

def filter_isec_data(isec_data, mode="genotype_differences"):
    """
    Filter isec data for visualization.
    """
    filtered_data = {}
    for chromosome, dfs in isec_data.items():
        filtered_dfs = []
        for df in dfs:
            if mode == "genotype_differences":
                filtered_df = df[df["Genotype"] != "."]
            elif mode == "missing_entries":
                filtered_df = df[df["Genotype"] == "."]
            else:
                logging.warning(f"Invalid mode: {mode}")
                filtered_df = pd.DataFrame()

            if not filtered_df.empty:
                filtered_dfs.append(filtered_df)

        if filtered_dfs:
            filtered_data[chromosome] = filtered_dfs

    return filtered_data

def create_combined_chromosome_page(chromosome_data, chromosome_lengths, samples, output_file):
    """
    Create a single-page visualization with a dark theme and vibrant colors.
    """
    logging.info("Creating combined visualization")
    try:
        # Validate output directory
        if not os.access(os.path.dirname(output_file), os.W_OK):
            logging.error(f"Output directory is not writable: {os.path.dirname(output_file)}")
            return

        # Basic setup
        num_chromosomes = len(chromosome_lengths)
        num_samples = len(samples)
        width = 2000
        height = max(200 * num_chromosomes + 300, 100 * num_samples + 300)
        margin = 150
        max_length = max(chromosome_lengths.values())
        x_scale = (width - 2 * margin) / max_length
        y_spacing = 200
        sample_offset = 10

        # Log basic dimensions and setup
        logging.info(f"Number of chromosomes: {num_chromosomes}")
        logging.info(f"Number of samples: {num_samples}")
        logging.info(f"Canvas dimensions: width={width}, height={height}")
        logging.info(f"x_scale: {x_scale}, y_spacing: {y_spacing}")

        # Ensure proper chromosome sorting (handles chromosome 01)
        sorted_chromosomes = sorted(chromosome_lengths.items(), key=lambda x: int(x[0]))

        # Initialize cairo surface
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        ctx = cairo.Context(surface)

        # Background
        ctx.rectangle(0, 0, width, height)
        ctx.set_source_rgb(0.1, 0.1, 0.1)
        ctx.fill()

        # Title
        ctx.set_source_rgb(1, 1, 1)
        ctx.set_font_size(28)
        ctx.move_to(margin, 50)
        ctx.show_text("Chromosome Visualization")

        # Colors for samples
        cmap = plt.cm.get_cmap("tab10")
        sample_colors = {sample: cmap(i / max(1, len(samples)))[:3] for i, sample in enumerate(samples)}
        logging.info(f"Sample colors assigned: {sample_colors}")

        # Legend
        legend_x = width - 300
        legend_y = 100
        ctx.set_source_rgb(1, 1, 1)
        ctx.set_font_size(24)
        ctx.move_to(legend_x, legend_y)
        ctx.show_text("Legend:")
        for i, (sample_name, color) in enumerate(sample_colors.items()):
            ctx.set_source_rgb(*color)
            ctx.rectangle(legend_x + 20, legend_y + (i + 1) * 30 - 10, 20, 20)
            ctx.fill()
            ctx.set_source_rgb(1, 1, 1)
            ctx.move_to(legend_x + 50, legend_y + (i + 1) * 30 + 5)
            ctx.show_text(sample_name)

        # Chromosome visualization
        y_position = 150
        ctx.move_to(margin - 50, 150)

        for idx, (chromosome, length) in enumerate(sorted_chromosomes):
            logging.info(f"Processing chromosome {chromosome} with length {length}")
            label_x = margin - 50
            label_y = y_position
            ctx.save()
            ctx.set_source_rgb(1, 1, 1)
            ctx.set_font_size(24)
            ctx.translate(label_x, label_y)
            ctx.rotate(-3.14159 / 2)
            ctx.show_text(f"{chromosome.zfill(2)}")  # Ensure chromosome labels like 01, 02, etc.
            ctx.restore()

            if chromosome in chromosome_data:
                samples_data = chromosome_data[chromosome]
                for sample_index, sample_name in enumerate(samples):
                    sample_data_list = [df for df in samples_data if df["Sample"].iloc[0] == sample_name]
                    if not sample_data_list:
                        logging.info(f"No data for sample {sample_name} on chromosome {chromosome}")
                        continue
                    sample_data = sample_data_list[0]
                    logging.info(f"Drawing data for sample {sample_name} on chromosome {chromosome}")

                    color = sample_colors[sample_name]
                    ctx.set_source_rgb(*color)
                    ctx.set_line_width(0.1)

                    sample_y_offset = y_position + (sample_index - num_samples / 2) * sample_offset
                    logging.info(f"Min: {sample_data['Position'].min()}, Max: {sample_data['Position'].max()}")
                    for _, row in sample_data.iterrows():
                        try:
                            pos = int(row["Position"])
                            x_pos = margin + pos * x_scale
                            ctx.arc(x_pos, sample_y_offset, 0.1, 0, 2 * 3.14159)
                            ctx.fill()
                        except Exception as e:
                            logging.warning(f"Error plotting position {row['Position']} for sample {sample_name}: {e}")

            y_position += y_spacing
            logging.info(f"Moving to next chromosome position: Y={y_position}")

        # Bottom axis with position bins
        ctx.set_source_rgb(1, 1, 1)
        ctx.set_font_size(20)
        bin_size = 200000
        for bin_start in range(0, max_length + bin_size, bin_size):
            x_pos = margin + bin_start * x_scale
            ctx.move_to(x_pos, height - 100)
            ctx.line_to(x_pos, height - 90)
            ctx.stroke()
            ctx.save()
            ctx.translate(x_pos, height - 70)
            ctx.rotate(-3.14159 / 2)  # Rotate 90 degrees for label
            ctx.show_text(f"{bin_start // 1000}K")
            ctx.restore()

        # Align all chromosomes to start uniformly
        ctx.set_source_rgb(1, 1, 1)
        ctx.set_line_width(2)
        ctx.move_to(margin, 150)
        ctx.line_to(margin, y_position)
        ctx.stroke()

        # Write output
        try:
            surface.write_to_png(output_file)
            logging.info(f"Visualization successfully saved to: {output_file}")
        except Exception as e:
            logging.error(f"Failed to write to file {output_file}: {e}")

    except Exception as e:
        logging.error(f"Error creating visualization: {e}")



def create_combined_visualizations_by_mode(isec_data, chromosome_lengths, samples, output_dir):
    """
    Create visualizations for genotype differences and missing entries.
    """
    genotype_diff_data = filter_isec_data(isec_data, mode="genotype_differences")
    genotype_diff_file = os.path.join(output_dir, "genotype_differences_visualization.png")
    create_combined_chromosome_page(genotype_diff_data, chromosome_lengths, samples, genotype_diff_file)

    missing_entries_data = filter_isec_data(isec_data, mode="missing_entries")
    missing_entries_file = os.path.join(output_dir, "missing_entries_visualization.png")
    create_combined_chromosome_page(missing_entries_data, chromosome_lengths, samples, missing_entries_file)

if __name__ == "__main__":
    isec_dir = "gvcf_comparison_results/isec"
    reference_genome_path = "qc/PlasmoDB-61_Pfalciparum3D7_Genome.fasta"
    output_dir = os.path.abspath("chromosome_visualizations")

    logging.info("Starting visualization script")

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        logging.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

    # Verify directory permissions
    if not os.access(output_dir, os.W_OK):
        logging.error(f"Output directory is not writable: {output_dir}")
        raise PermissionError(f"Cannot write to output directory: {output_dir}")

    chromosome_lengths = get_chromosome_lengths(reference_genome_path)
    isec_data, samples = parse_isec_vcfs(isec_dir)

    create_combined_visualizations_by_mode(isec_data, chromosome_lengths, samples, output_dir)

    logging.info("Script completed successfully")

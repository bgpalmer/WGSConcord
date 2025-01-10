#!/bin/bash

# Directories
TERRA_DIR="Terra/gcvf_downloads"
KARAMOKO_DIR="Karamoko"
OUTPUT_DIR="gvcf_comparison_results"
STATS_DIR="$OUTPUT_DIR/stats"
VARIANTS_DIR="$OUTPUT_DIR/variants"
ISEC_DIR="$OUTPUT_DIR/isec"
GENOME="qc/PlasmoDB-61_Pfalciparum3D7_Genome.fasta"

# Create necessary directories
mkdir -p "$OUTPUT_DIR" "$STATS_DIR" "$VARIANTS_DIR" "$ISEC_DIR"

# Extract sample name from Terra files (up to `.haplotype_caller`)
extract_terra_sample_key() {
    basename "$1" | sed -E 's/^(.*)\.haplotype_caller.*/\1/'
}

# Map Terra chromosome names to Karamoko (remove zero-fill)
map_chromosome_name() {
    local chrom_name="$1"
    echo "${chrom_name/Pf3D7_/chr}" | sed 's/_v3//' | sed 's/^chr0/chr/'
}

# Clean VCF files by removing `<NON_REF>` lines
clean_vcf_file() {
    local input_vcf="$1"
    local output_vcf="$2"

    echo "Cleaning VCF file: $input_vcf -> $output_vcf"
    bcftools sort "$input_vcf" -Oz -o "sorted_${input_vcf##*/}"
    bcftools norm -m - "sorted_${input_vcf##*/}" | egrep -v "<NON_REF>" | bcftools convert -Oz -o "$output_vcf"
    bcftools index -f "$output_vcf"
    echo "Cleaned VCF file written to: $output_vcf"
}

# Split Terra GVCFs by chromosome
split_terra_by_chromosome() {
    local terra_file="$1"
    local output_dir="$2"
    mkdir -p "$output_dir"
    for chrom in $(seq 1 14); do
        chrom_name=$(printf "Pf3D7_%02d_v3" "$chrom")
        output_file="$output_dir/${chrom_name}_snps_indels.g.vcf.gz"

        if [ ! -f "$output_file" ]; then
            echo "Splitting Terra file $terra_file for chromosome $chrom_name..."
            if [ ! -f "${terra_file}.tbi" ]; then
                echo "Index missing for $terra_file. Generating index..."
                bcftools index -f "$terra_file"
            fi

            bcftools view -v snps -r "$chrom_name" -Oz -o "$output_file" "$terra_file"
            bcftools index -f "$output_file"
        fi
    done
}

# Search for Karamoko file matching the sample and chromosome
find_karamoko_file() {
    local sample_key="$1"
    local chrom_name="$2"
    find "$KARAMOKO_DIR" -type f -name "*${sample_key}*${chrom_name}.g.vcf.gz" | head -n 1
}

# Main Processing Loop
for TERRA_FILE in $(find "$TERRA_DIR" -type f -name "*.g.vcf.gz" | sort); do
    TERRA_SAMPLE_KEY=$(extract_terra_sample_key "$TERRA_FILE")
    if [ -z "$TERRA_SAMPLE_KEY" ]; then
        echo "Error: Could not extract sample key from $TERRA_FILE. Skipping."
        continue
    fi

    echo "Processing Terra sample: $TERRA_SAMPLE_KEY"

    # Create a sample-specific stats directory
    SAMPLE_STATS_DIR="$STATS_DIR/$TERRA_SAMPLE_KEY"
    mkdir -p "$SAMPLE_STATS_DIR"

    # Split Terra file by chromosome
    TERRA_CHROM_DIR="$VARIANTS_DIR/${TERRA_SAMPLE_KEY}_terra_chromosomes"
    split_terra_by_chromosome "$TERRA_FILE" "$TERRA_CHROM_DIR"

    # Process each chromosome
    for chrom in $(seq 1 14); do
        chrom_name=$(printf "Pf3D7_%02d_v3" "$chrom")
        karamoko_chrom_name=$(map_chromosome_name "$chrom_name")
        terra_chrom_file="$TERRA_CHROM_DIR/${chrom_name}_snps_indels.g.vcf.gz"

        # Clean Terra VCF
        cleaned_terra_vcf="${terra_chrom_file%.vcf.gz}_cleaned.vcf.gz"
        clean_vcf_file "$terra_chrom_file" "$cleaned_terra_vcf"

        # Find corresponding Karamoko file
        karamoko_file=$(find_karamoko_file "$TERRA_SAMPLE_KEY" "$karamoko_chrom_name")
        if [ -z "$karamoko_file" ]; then
            echo "No matching Karamoko file found for $TERRA_SAMPLE_KEY and chromosome $karamoko_chrom_name. Skipping."
            continue
        fi

        echo "Found matching Karamoko file: $karamoko_file"

        # Clean Karamoko VCF
        cleaned_karamoko_vcf="${karamoko_file%.vcf.gz}_cleaned.vcf.gz"
        clean_vcf_file "$karamoko_file" "$cleaned_karamoko_vcf"

        # Perform comparisons with bcftools isec
        isec_output_dir="${ISEC_DIR}/${TERRA_SAMPLE_KEY}_${chrom_name}"
        mkdir -p "$isec_output_dir"
        echo "Running bcftools isec..."
        bcftools isec -p "$isec_output_dir" "$cleaned_terra_vcf" "$cleaned_karamoko_vcf"

        # Generate stats for original and comparison files
        echo "Generating bcftools stats..."
        bcftools stats "$cleaned_terra_vcf" > "$SAMPLE_STATS_DIR/${chrom_name}_terra.stats.txt"
        bcftools stats "$cleaned_karamoko_vcf" > "$SAMPLE_STATS_DIR/${chrom_name}_karamoko.stats.txt"

        for isec_file in "$isec_output_dir"/*.vcf; do
            isec_stats_file="$SAMPLE_STATS_DIR/$(basename "${isec_file%.vcf}")_stats.txt"
            bcftools stats "$isec_file" > "$isec_stats_file"
            echo "Generated stats for $isec_file -> $isec_stats_file"
        done
    done
done

echo "All comparisons and stats generation completed. Outputs are in $OUTPUT_DIR."

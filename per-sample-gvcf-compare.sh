#!/bin/bash

# Directories
TERRA_DIR="Terra/gcvf_downloads"
KARAMOKO_DIR="Karamoko"
OUTPUT_DIR="gvcf_comparison_results"
STATS_DIR="$OUTPUT_DIR/stats"
PLOTS_DIR="$OUTPUT_DIR/plots"
VARIANTS_DIR="$OUTPUT_DIR/variants"

mkdir -p "$OUTPUT_DIR" "$STATS_DIR" "$PLOTS_DIR" "$VARIANTS_DIR"

# Function to extract identifier before `.haplotype`
extract_sample_key() {
    basename "$1" | sed -E 's/^(.*)\.haplotype.*/\1/'
}

# Function to generate index files for Terra GVCFs
generate_terra_indices() {
    local terra_file="$1"
    echo "Checking index for ${terra_file}..."
    if [ ! -f "${terra_file}.tbi" ]; then
        echo "Index file missing for ${terra_file}. Generating..."
        if bcftools index -f "$terra_file"; then
            echo "Index file generated successfully for ${terra_file}."
        else
            echo "Failed to generate index for ${terra_file}. Skipping this file."
            return 1
        fi
    fi
    return 0
}

# Function to split Terra GVCF by chromosome
split_terra_by_chromosome() {
    local terra_file="$1"
    local output_dir="$2"
    mkdir -p "$output_dir"
    for chrom in $(seq 1 14); do
        chrom_name="Pf3D7_${chrom}_v3"
        output_file="$output_dir/${chrom_name}.vcf.gz"

        if [ ! -f "$output_file" ]; then
            echo "Extracting chromosome ${chrom_name} from $terra_file..."
            if [ ! -f "${terra_file}.tbi" ]; then
                echo "Missing index file for $terra_file. Attempting to regenerate..."
                generate_terra_indices "$terra_file" || return 1
            fi

            if bcftools view -r "${chrom_name}" -Oz -o "$output_file" "$terra_file"; then
                bcftools index -f "$output_file"
                echo "Chromosome ${chrom_name} extracted successfully."
            else
                echo "Error extracting chromosome ${chrom_name} from $terra_file. Skipping this chromosome."
                continue
            fi
        fi
    done
    return 0
}

# Function to map Terra chromosome names to Karamoko names
map_chromosome_name() {
    local chrom_name="$1"
    echo "${chrom_name/Pf3D7_/chr}" | sed 's/_v3//' | sed 's/^chr0/chr/'  # Remove zero-fill
}

# Function to compress and index GVCFs
ensure_bgzip_and_index() {
    local file="$1"
    # if [[ "${file}" != *.gz ]]; then
    #     echo "Compressing ${file} with BGZIP..."
    #     if bgzip -f "$file"; then
    #         echo "Compressed ${file} successfully."
    #     else
    #         echo "Failed to compress ${file}. Skipping."
    #         return 1
    #     fi
    # fi

    local gz_file="${file}"
    if [ ! -f "${gz_file}.tbi" ]; then
        echo "Generating index for ${gz_file}..."
        if bcftools index -f "$gz_file"; then
            echo "Index file generated successfully for ${gz_file}."
        else
            echo "Failed to generate index for ${gz_file}. Skipping."
            return 1
        fi
    fi
    return 0
}

# Function to search for Karamoko GVCFs
find_karamoko_files() {
    local sample_key="$1"
    local chrom_name="$2"

    # Dynamically find batch subdirectories
    for batch_dir in "$KARAMOKO_DIR"/*_GVCF_results; do
        local chrom_dir="$batch_dir/$chrom_name"
        if [ -d "$chrom_dir" ]; then
            find "$chrom_dir" -type f -name "*${sample_key}*.g.vcf.gz" | sort
        fi
    done
}

# Main Processing Loop
for TERRA_FILE in $(find "$TERRA_DIR" -type f -name "*.g.vcf.gz" | sort); do
    SAMPLE_KEY=$(extract_sample_key "$TERRA_FILE")
    if [ -z "$SAMPLE_KEY" ]; then
        echo "Could not extract sample key from $TERRA_FILE. Skipping."
        continue
    fi

    echo "Processing sample: $SAMPLE_KEY"

    # Generate index if missing
    generate_terra_indices "$TERRA_FILE" || continue

    # Split Terra GVCF by chromosome
    TERRA_CHROM_DIR="$VARIANTS_DIR/${SAMPLE_KEY}_terra_chromosomes"
    split_terra_by_chromosome "$TERRA_FILE" "$TERRA_CHROM_DIR" || continue

    # Compare each chromosome
    for chrom in $(seq 1 14); do
        chrom_name="Pf3D7_${chrom}_v3"
        terra_variants="$TERRA_CHROM_DIR/${chrom_name}.vcf.gz"
        karamoko_chrom_name=$(map_chromosome_name "$chrom_name")

        echo "Searching for Karamoko GVCFs for chromosome $karamoko_chrom_name..."
        echo "SAMPLE_KEY: $SAMPLE_KEY"
        karamoko_variants=$(find_karamoko_files "$SAMPLE_KEY" "$karamoko_chrom_name")

        echo $karamoko_variants

        if [ -z "$karamoko_variants" ]; then
            echo "No Karamoko GVCF found for chromosome ${karamoko_chrom_name}. Skipping."
            continue
        fi

        # Ensure Karamoko GVCFs are compressed and indexed
        ensure_bgzip_and_index "$karamoko_variants" || continue

        # Perform statistics comparison
        sample_stats_file="$STATS_DIR/${SAMPLE_KEY}_${chrom_name}.stats.txt"
        echo "Generating statistics for sample $SAMPLE_KEY chromosome ${chrom_name}..."
        if bcftools stats --threads 4 "$terra_variants" "${karamoko_variants}" > "$sample_stats_file"; then
            echo "Statistics saved for $SAMPLE_KEY chromosome ${chrom_name}."
        else
            echo "Error generating statistics for $SAMPLE_KEY chromosome ${chrom_name}."
            continue
        fi
    done

    # Generate precision and recall plot
    python3 generate_precision_recall_plot.py "$SAMPLE_KEY" "$STATS_DIR" "$PLOTS_DIR/${SAMPLE_KEY}_precision_recall.png"
done

echo "All comparisons completed. Outputs are in $OUTPUT_DIR."

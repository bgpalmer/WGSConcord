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

# Function to format chromosome names with zero padding for 1-9
format_chromosome_name() {
    local chrom_num="$1"
    printf "Pf3D7_%02d_v3" "$chrom_num"  # Adds zero padding for numbers 1-9
}

# Function to split and filter Terra GVCF by chromosome
split_and_filter_terra_by_chromosome() {
    local terra_file="$1"
    local output_dir="$2"
    mkdir -p "$output_dir"
    for chrom in $(seq 1 14); do
        chrom_name=$(format_chromosome_name "$chrom")  # Use the formatted chromosome name
        output_file="$output_dir/${chrom_name}_snps_indels.g.vcf.gz"

        if [ ! -f "$output_file" ]; then
            echo "Extracting and filtering SNPs and indels for chromosome ${chrom_name} from $terra_file..."
            if [ ! -f "${terra_file}.tbi" ]; then
                echo "Missing index file for $terra_file. Attempting to regenerate..."
                generate_terra_indices "$terra_file" || return 1
            fi

            if bcftools view -v snps,indels -r "${chrom_name}" -Oz -o "$output_file" "$terra_file"; then
                bcftools index -f "$output_file"
                echo "Chromosome ${chrom_name} SNPs and indels extracted successfully."
            else
                echo "Error extracting SNPs and indels for chromosome ${chrom_name} from $terra_file. Skipping this chromosome."
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

# Function to filter Karamoko GVCF by SNPs and indels
filter_karamoko_by_chromosome() {
    local karamoko_file="$1"
    local output_file="$2"

    if [ ! -f "$karamoko_file" ]; then
        echo "Error: Karamoko file $karamoko_file not found. Skipping."
        return 1
    fi

    if [ ! -f "$output_file" ]; then
        echo "Filtering SNPs and indels for $karamoko_file..."
        if [ ! -f "${karamoko_file}.tbi" ]; then
            echo "Missing index file for $karamoko_file. Generating index..."
            ensure_bgzip_and_index "$karamoko_file" || return 1
        fi

        if bcftools view -v snps,indels -Oz -o "$output_file" "$karamoko_file"; then
            bcftools index -f "$output_file"
            echo "Filtered SNPs and indels written to $output_file."
        else
            echo "Error filtering SNPs and indels for $karamoko_file. Skipping."
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

# Function to compress and index GVCFs
ensure_bgzip_and_index() {
    local file="$1"

    # # Compress with bgzip if not already compressed
    # if [[ "${file}" != *.gz ]]; then
    #     echo "Compressing ${file} with BGZIP..."
    #     if bgzip -f "$file"; then
    #         echo "Compressed ${file} successfully."
    #     else
    #         echo "Failed to compress ${file}. Skipping."
    #         return 1
    #     fi
    # fi

    local gz_file="${file%.gz}.gz"
    # Index the gzipped file
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


# Main Processing Loop
for TERRA_FILE in $(find "$TERRA_DIR" -type f -name "*.g.vcf.gz" | sort); do
    SAMPLE_KEY=$(extract_sample_key "$TERRA_FILE")
    if [ -z "$SAMPLE_KEY" ]; then
        echo "Error: Could not extract sample key from $TERRA_FILE. Skipping."
        continue
    fi

    echo "Processing sample: $SAMPLE_KEY"

    # Generate index for Terra file if missing
    generate_terra_indices "$TERRA_FILE" || continue

    # Split Terra file by chromosome
    TERRA_CHROM_DIR="$VARIANTS_DIR/${SAMPLE_KEY}_terra_chromosomes"
    split_and_filter_terra_by_chromosome "$TERRA_FILE" "$TERRA_CHROM_DIR" || continue

    # Process each chromosome
    for chrom in $(seq 1 14); do
        chrom_name=$(format_chromosome_name "$chrom")
        terra_variants="$TERRA_CHROM_DIR/${chrom_name}_snps_indels.g.vcf.gz"
        karamoko_chrom_name=$(map_chromosome_name "$chrom_name")

        # Find corresponding Karamoko file
        karamoko_variants=$(find_karamoko_files "$SAMPLE_KEY" "$karamoko_chrom_name")
        if [ -z "$karamoko_variants" ]; then
            echo "No Karamoko file found for chromosome $karamoko_chrom_name. Skipping."
            continue
        fi

        # Filter Karamoko file
        KARAMOKO_CHROM_DIR="$VARIANTS_DIR/${SAMPLE_KEY}_karamoko_chromosomes"
        mkdir -p "$KARAMOKO_CHROM_DIR"
        karamoko_filtered="$KARAMOKO_CHROM_DIR/${karamoko_chrom_name}_snps_indels.g.vcf.gz"

        filter_karamoko_by_chromosome "$karamoko_variants" "$karamoko_filtered" || continue

        # Generate intersection/complement sets using bcftools isec
        INTERSECT_DIR="$OUTPUT_DIR/intersect/${SAMPLE_KEY}_${chrom_name}"
        mkdir -p "$INTERSECT_DIR"
        echo "Creating intersection/complement sets for $chrom_name..."
        bcftools isec -p "$INTERSECT_DIR" -Oz "$terra_variants" "$karamoko_filtered"

        # Run bcftools stats on Terra and Karamoko files
        echo "Generating bcftools stats for $chrom_name..."
        terra_stats_file="$STATS_DIR/${SAMPLE_KEY}_${chrom_name}_terra.stats.txt"
        karamoko_stats_file="$STATS_DIR/${SAMPLE_KEY}_${chrom_name}_karamoko.stats.txt"

        bcftools stats "$terra_variants" > "$terra_stats_file"
        bcftools stats "$karamoko_filtered" > "$karamoko_stats_file"

        # Generate bcftools stats on intersection sets
        intersect_stats_file="$STATS_DIR/${SAMPLE_KEY}_${chrom_name}_intersect.stats.txt"
        if bcftools stats "$INTERSECT_DIR/0003.vcf.gz" > "$intersect_stats_file"; then
            echo "Intersection statistics generated for $chrom_name."
        else
            echo "Error generating intersection statistics for $chrom_name."
            continue
        fi

        # Generate precision and recall plot
        echo "Generating precision and recall plot..."
        python3 generate_precision_recall_plot.py \
            "$SAMPLE_KEY" "$STATS_DIR" "$PLOTS_DIR/${SAMPLE_KEY}_${chrom_name}_precision_recall.png"
    done
done

echo "All comparisons completed. Outputs are in $OUTPUT_DIR."

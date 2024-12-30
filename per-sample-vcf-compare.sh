#!/bin/bash

# Define directories
TERRA_VCF="Terra/submissions_1717f5fa-7b6c-4312-95b0-801562a1bb6b_SRWholeGenome_0fae6438-5610-4b17-8a39-3b4ec2552d28_call-ScoreIndelVariantAnnotations_S3_WGS-batchC_4051002301-ST140_S9_ALL_scored.vcf"
KARAMOKO_DIR="Karamoko/20230808_batchC_GVCF_results"

# Output directories
OUTPUT_DIR="vcf_comparison_results"
STATS_DIR="$OUTPUT_DIR/stats"
PLOTS_DIR="$OUTPUT_DIR/plots"
ISEC_DIR="$OUTPUT_DIR/isec"
PDF_REPORT="$OUTPUT_DIR/report.pdf"
COMBINED_VCF="$OUTPUT_DIR/combined_karamoko.vcf.gz"
CONCORDANCE_OUTPUT="$OUTPUT_DIR/concordance_summary.txt"

mkdir -p "$OUTPUT_DIR" "$STATS_DIR" "$PLOTS_DIR" "$ISEC_DIR"

# Compress and index the Terra VCF if not already compressed
if [[ "$TERRA_VCF" != *.gz ]]; then
    echo "Compressing Terra VCF..."
    bgzip -c "$TERRA_VCF" > "${TERRA_VCF}.gz" || { echo "Error compressing $TERRA_VCF"; exit 1; }
    TERRA_VCF="${TERRA_VCF}.gz"
fi

# Index the Terra VCF if index doesn't exist
if [ ! -f "${TERRA_VCF}.csi" ]; then
    echo "Indexing Terra VCF..."
    bcftools index "$TERRA_VCF" || { echo "Error indexing $TERRA_VCF"; exit 1; }
fi

# Combine Karamoko VCF files into one
echo "Combining Karamoko VCF files into a single file..."
if [ ! -f "$COMBINED_VCF" ]; then
    bcftools concat -Oz -o "$COMBINED_VCF" "$KARAMOKO_DIR"/chr*/S3_WGS-batchC_4051002301-ST140_S9.*.g.vcf \
        && bcftools index "$COMBINED_VCF" || { echo "Error concatenating or indexing Karamoko VCF files"; exit 1; }
else
    echo "Combined VCF already exists: $COMBINED_VCF"
fi

# Perform intersection analysis using bcftools isec
echo "Performing intersection analysis using bcftools isec..."
bcftools isec -p "$ISEC_DIR" -n=2 "$TERRA_VCF" "$COMBINED_VCF" \
    && echo "Intersection analysis results saved to $ISEC_DIR" \
    || { echo "Error performing intersection analysis"; exit 1; }

# Summarize intersection results
SHARED_COUNT=$(grep -vc '^#' "$ISEC_DIR/0002.vcf")
TERRA_UNIQUE_COUNT=$(grep -vc '^#' "$ISEC_DIR/0000.vcf")
KARAMOKO_UNIQUE_COUNT=$(grep -vc '^#' "$ISEC_DIR/0001.vcf")
echo -e "Intersection Summary:\nShared Variants: $SHARED_COUNT\nUnique to Terra: $TERRA_UNIQUE_COUNT\nUnique to Karamoko: $KARAMOKO_UNIQUE_COUNT" > "$OUTPUT_DIR/intersection_summary.txt"

# Perform comparison with bcftools stats
echo "Generating statistics for Terra VCF and combined Karamoko VCF..."
STATS_FILE="$STATS_DIR/combined_stats.txt"
bcftools stats -s - "$TERRA_VCF" "$COMBINED_VCF" > "$STATS_FILE" \
    && echo "Statistics saved to $STATS_FILE" || { echo "Error generating statistics"; exit 1; }

# Extract all concordance data
echo "Extracting concordance data..."
CONCORDANCE_SECTIONS=("GCsAF" "NRDs" "GCiAF" "NRDi" "GCsS" "GCiS" "GCTs" "GCTi")
: > "$CONCORDANCE_OUTPUT" # Clear existing file or create a new one
for SECTION in "${CONCORDANCE_SECTIONS[@]}"; do
    echo "### $SECTION ###" >> "$CONCORDANCE_OUTPUT"
    grep "^$SECTION" "$STATS_FILE" >> "$CONCORDANCE_OUTPUT"
    echo >> "$CONCORDANCE_OUTPUT"
done
echo "All concordance data saved to $CONCORDANCE_OUTPUT"

# Generate plots from statistics
echo "Generating plots from statistics..."
plot-vcfstats "$STATS_FILE" -p "$PLOTS_DIR" --no-PDF \
    && echo "Plots saved to $PLOTS_DIR" || { echo "Error generating plots"; exit 1; }

# Create a consolidated PDF report from plots and concordance summary
echo "Generating PDF report..."
SUMMARY_TEXT_PDF="$PLOTS_DIR/summary_text.pdf"
if command -v pandoc &>/dev/null; then
    pandoc "$CONCORDANCE_OUTPUT" -o "$SUMMARY_TEXT_PDF" \
        && echo "Concordance summary converted to PDF" \
        || echo "Error generating summary text PDF. Ensure Pandoc is installed."
else
    echo "Pandoc is required to convert text summary to PDF."
fi

if command -v convert &>/dev/null; then
    convert "$PLOTS_DIR"/*.png "$SUMMARY_TEXT_PDF" "$PDF_REPORT" \
        && echo "PDF report saved to $PDF_REPORT" \
        || echo "Error generating PDF report. Ensure ImageMagick is installed."
else
    echo "ImageMagick (convert) is required to generate the PDF report."
fi

echo "All processes completed successfully."
echo "Concordance summary available at: $CONCORDANCE_OUTPUT"
echo "Intersection summary available at: $OUTPUT_DIR/intersection_summary.txt"

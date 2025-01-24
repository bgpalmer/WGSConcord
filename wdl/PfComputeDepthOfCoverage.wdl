version 1.0

workflow PfReadDepthWorkflow {
    input {
        String sample_id
        File pf_bam
        File pf_bam_index
        File ref_map_file  # Accepts the ref_map file as input
        File regions_bed_file  # BED file specifying regions of interest
        Int task_memory_gb = 16  # Set default memory allocation (adjustable)
        String output_dir = "output"  # Default output directory
    }

    # Read the reference map file into a Map[String, String]
    Map[String, String] ref_map = read_map(ref_map_file)

    call CalculateReadDepth {
        input:
            sample_id = sample_id,
            pf_bam = pf_bam,
            pf_bam_index = pf_bam_index,
            ref_fasta = ref_map["fasta"],
            ref_dict = ref_map["dict"],
            ref_fasta_index = ref_map["fai"],
            regions_bed = regions_bed_file,
            task_memory_gb = task_memory_gb
    }

    output {
        File read_coverage_tsv = CalculateReadDepth.output_tsv
    }
}

task CalculateReadDepth {
    input {
        String sample_id
        File pf_bam
        File pf_bam_index
        File ref_fasta
        File ref_dict
        File ref_fasta_index
        File regions_bed  # New input for regions BED file
        Int task_memory_gb
    }

    command <<<
        mkdir -p TMP

        gatk --java-options "-Xmx~{task_memory_gb}g -Xms~{task_memory_gb}g" DepthOfCoverage \
            -R ~{ref_fasta} \
            -I ~{pf_bam} \
            -L ~{regions_bed} \
            --output-format TABLE \
            --omit-locus-table true \
            -O depth_of_coverage \
            --tmp-dir TMP

        # Post-process the output to ensure formatting matches expectations
        awk -F"\t" -v OFS="\t" '{ print $0 }' depth_of_coverage.sample_summary > tmp.sample_summary
        sed '${/^Total/d;}' tmp.sample_summary > filtered.sample_summary
        awk -F"\t" -v OFS="\t" '{ print $0 }' filtered.sample_summary > ~{sample_id}_read_coverage.tsv
    >>>

    output {
        File output_tsv = "~{sample_id}_read_coverage.tsv"
    }

    runtime {
        memory: "~{task_memory_gb}G"
        disks: "local-disk 250 SSD"
        cpu: 1
        docker: "broadinstitute/gatk:4.6.1.0"
    }
}

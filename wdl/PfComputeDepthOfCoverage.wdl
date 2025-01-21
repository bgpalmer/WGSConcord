version 1.0

workflow PfReadDepthWorkflow {
    input {
        String sample_id
        File pf_bam
        File pf_bam_index
        File ref_map_file  # Accepts the ref_map file as input
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
        Int task_memory_gb
    }

    command <<<
        for i in {01..14}; do
            gatk --java-options "-Xmx~{task_memory_gb}g -Xms~{task_memory_gb}g" DepthOfCoverage \
                -R ~{ref_fasta} \
                -O chr${i} \
                --output-format TABLE \
                -L Pf3D7_${i}_v3 \
                --omit-locus-table true \
                -I ~{pf_bam} \
            sed '${/^Total/d;}' chr${i}.sample_summary > tmp.chr${i}.sample_summary
            awk -F"\t" -v OFS="\t" '{ print $0, $(NF) = "chr'${i}'" }' tmp.chr${i}.sample_summary > chr${i}.sample2_summary
        done

        cat *.sample2_summary | awk '!/sample_id/ {print $0}' | sed '1isample_id\ttotal\tmean\tthird_quartile\tmedian\tfirst_quartile\tbases_perc_above_15\tchromosome' > ~{sample_id}_read_coverage.tsv
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

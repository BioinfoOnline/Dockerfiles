#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//process1: fastqc
process quality_control {
    container 'Dockerfile'      //one Dockerfile for all the required tools

    input:
    path reads

    output:
    path("${params.outputDir}/qc/") into qc_results

    script:
    """
    mkdir -p ${params.outputDir}/qc/
    fastqc ${reads} -o ${params.outputDir}/qc/
    """
}

//process 2: trimming of adaptors
process trim_adapters {
    container 'Dockerfile'

    input:
    path reads from quality_control.out

    output:
    path("${params.outputDir}/trimmed/") into trimmed_results

    script:
    """
    mkdir -p ${params.outputDir}/trimmed/
    java -jar /opt/trimmomatic-0.39.jar PE -phred33 ${reads} ${reads} \
        ${params.outputDir}/trimmed/trimmed_paired.fq.gz ${params.outputDir}/trimmed/trimmed_unpaired.fq.gz \
        ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

//process3: indexing the ref
process index_reference {
    container 'Dockerfile'

    input:
    path reference

    output:
    path("${params.outputDir}/reference/") into indexed_results

    script:
    """
    mkdir -p ${params.outputDir}/reference/
    bwa index -p ${params.outputDir}/reference/index ${reference}
    """
}

//process4: alignment
process align_reads {
    container 'Dockerfile'

    input:
    path trimmed_reads from trim_adapters.out
    path index_prefix from index_reference.out

    output:
    path("${params.outputDir}/aligned/") into alignment_results

    script:
    """
    mkdir -p ${params.outputDir}/aligned/
    bwa mem ${index_prefix}/index ${trimmed_reads} > ${params.outputDir}/aligned/aligned.sam
    """
}

//process 5 deduplication
process deduplicate {
    container 'Dockerfile'

    input:
    path aligned_sam from align_reads.out

    output:
    path("${params.outputDir}/dedup/") into deduped_results

    script:
    """
    mkdir -p ${params.outputDir}/dedup/

    samtools view -Sb ${aligned_sam} > ${params.outputDir}/dedup/aligned.bam    # Convert SAM to BAM

    samtools markdup -r ${params.outputDir}/dedup/aligned.bam ${params.outputDir}/dedup/aligned_dedup.bam       # Mark duplicates

    samtools index ${params.outputDir}/dedup/aligned_dedup.bam      # Index the BAM file
    """
}

//process 6: variant call
process variant_calling {
    container 'Dockerfile'

    input:
    path aligned_bam
    path reference
    path output_dir

    output:
    path("${output_dir}/variant_calls/") into variant_calls

    script:
    """
    mkdir -p ${output_dir}/variant_calls/

    # Run GATK HaplotypeCaller for variant calling
    gatk HaplotypeCaller -R ${reference} -I ${aligned_bam} -O ${output_dir}/variant_calls/raw_variants.vcf
    """
}

// Workflow
workflow {
    // Input parameters
    params.reads = file(params.reads ?: 'data/*.fastq.gz')
    params.reference = file(params.reference ?: 'data/sequence.fasta')  // Reference file is 'sequence.fasta'
    params.outputDir = params.outputDir ?: 'results'

    // Ensure output directory exists
    process.mkdir = true

    // Processes
    qc = quality_control(params.reads, params.outputDir)
    trimmed = trim_adapters(qc.out, params.outputDir)
    indexed = index_reference(params.reference, params.outputDir)
    aligned = align_reads(trimmed.out, indexed.out, params.outputDir)
    deduped = deduplicate(aligned.out, params.outputDir)
    variant_calling(deduped.out, params.reference, params.outputDir)
}

// in terminal
// making changes in the .config file for making it run faster
//nextflow run variant_calling.nf -with-docker
//nextflow run variant_calling.nf --reads /home/vidhi/project_4th_sem/Docker_var_calling/data/*.fastq.gz --reference /home/vidhi/project_4th_sem/Docker_var_calling/data/sequence.fasta --outputDir /home/vidhi/project_4th_sem/Docker_var_calling/data/results

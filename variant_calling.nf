#!/usr/bin/env nextflow

nextflow.enable.dsl=2


workflow {
    params.reads = file(params.reads ?: 'data/*.fastq')  
    params.reference = file(params.reference ?: 'data/*.fasta') 
    params.outputDir = params.outputDir ?: 'results'
   

    qc_result = quality_control(params.reads)
    read1 = file("/home/vidhi/project_4th_sem/check/data/SRR12345678_1.fastq")
    read2 = file("/home/vidhi/project_4th_sem/check/data/SRR12345678_2.fastq")
    //output_dir = file("/home/vidhi/project_4th_sem/check/results/trimmed/")

    trimmed_result = trim_adapters(read1, read2)
    //trimmed_result = trim_adapters(qc_result)  
    index_result = index_reference(params.reference)  
    aligned_result = align_reads(trimmed_result, index_result)  
    deduped_result = deduplicate(aligned_result)  
    variant_calling_result = variant_calling(deduped_result, params.reference)  
    variant_annotation(variant_calling_result, params.reference)  
}

//process1: fastqc
process quality_control {
    container 'check_image' 

    input:
    path reads

    output:
    path("${params.outputDir}/qc/*_fastqc.zip")

    script:
    """
    mkdir -p ${params.outputDir}/qc/
    echo "Directory exists"
    fastqc ${reads} -o ${params.outputDir}/qc/
    """
}

//process2: trimming
process trim_adapters {
    container 'check_image'

    input:
    path read1
    path read2
    
    
    output:
    path("${params.outputDir}/trimmed/")

    script:
    """
    mkdir -p ${params.outputDir}/trimmed/
    echo "Directory exists."

    BASE_NAME=\$(basename ${read1} .fastq)
    
    java -jar /usr/local/bin/trimmomatic.jar PE -phred33 \
        ${read1} ${read2} \
        ${params.outputDir}/\${BASE_NAME}_R1_paired.fastq ${params.outputDir}/\${BASE_NAME}_R1_unpaired.fastq \
        ${params.outputDir}/\${BASE_NAME}_R2_paired.fastq ${params.outputDir}/\${BASE_NAME}_R2_unpaired.fastq \
        ILLUMINACLIP:/opt/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    echo "Adaptor Trimming is done."
    """
}

//indexing the reference
process index_reference {
    container 'check_image'

    input:
    path reference

    output:
    path("${params.outputDir}/reference/")

    script:
    """
    mkdir -p ${params.outputDir}/reference/
    echo "Directory ${params.outputDir}/reference/ exists"
    
    bowtie2-build ${reference} ${params.outputDir}/reference/reference_index
    echo "Reference is Indexed."
    """
}

//alignment
process align_reads {
    container 'check_image'

    input:
    path trimmed_result
    path index_result

    output:
    path("${params.outputDir}/aligned/")

    script:
    """
    mkdir -p ${params.outputDir}/aligned/
    echo "Directory ${params.outputDir}/aligned/ exists"
    
    # need to specify the paired reads for trimming
    trimmed1=\$(echo ${trimmed_result} | cut -d' ' -f1)
    trimmed2=\$(echo ${trimmed_result} | cut -d' ' -f2)

    bowtie2 -x ${index_result}/reference_index -1 \${trimmed1} -2 \${trimmed2} -S ${params.outputDir}/aligned/aligned.sam
    echo "Alignment is done."
    """
}

// deduplication
process deduplicate {
    container 'check_image'
 
    input:
    path aligned_sam

    output:
    path("${params.outputDir}/dedup/")

    script:
    """
    mkdir -p ${params.outputDir}/dedup/
    echo "Directory ${params.outputDir}/dedup/ exists"
    samtools view -Sb ${aligned_sam}/*.sam > ${params.outputDir}/dedup/aligned.bam		#sam to bam
    samtools markdup -r ${params.outputDir}/dedup/aligned.bam ${params.outputDir}/dedup/aligned_dedup.bam		#deduplicating
    samtools index ${params.outputDir}/dedup/aligned_dedup.bam		#indexing the unique reads
    """
}

// var. calling
process variant_calling {
    container 'check_image'

    input:
    path aligned_bam
    path reference

    output:
    path("${params.outputDir}/variant_calls/")

    script:
    """
    mkdir -p ${params.outputDir}/variant_calls/
    echo "Directory ${params.outputDir}/variant_calls/ exists"
    
    ls -l ${aligned_bam}
    ls -l ${reference}
    
    samtools mpileup -f ${reference} ${aligned_bam} > ${params.outputDir}/variant_calls/pileup.txt
    echo "Pileup generated."

    # Call SNPs and indels using VarScan
    java -jar /usr/local/bin/varscan.jar mpileup2snp ${params.outputDir}/variant_calls/pileup.txt --output-vcf 1 > ${params.outputDir}/variant_calls/raw_variants.vcf
    echo "SNPs called."

    # Call Indels using VarScan
    java -jar /usr/local/bin/varscan.jar mpileup2indel ${params.outputDir}/variant_calls/pileup.txt --output-vcf 1 >> ${params.outputDir}/variant_calls/raw_variants.vcf
    echo "Indels called."
    """

}

//annotation
process variant_annotation {
    container 'check_image'

    input:
    path vcf_file
    path reference

    output:
    path("${oparams.outputDir}/annotated_variants/")

    script:
    """
    mkdir -p ${params.outputDir}/annotated_variants/
    echo "Directory ${params.outputDir}/annotated_variants exists"
    
    snpEff -v -c /opt/snpEff/snpEff.config -genome ${reference} ${vcf_file}/*.vcf > ${params.outputDir}/annotated_variants/annotated_variants.vcf
    echo "Annotation is done"
    """
}



// in terminal
// making changes in the .config file for making it run faster
//nextflow run variant_calling.nf -with-docker
//nextflow run variant_calling.nf --reads /home/vidhi/project_4th_sem/Docker_var_calling/data/*.fastq.gz --reference /home/vidhi/project_4th_sem/Docker_var_calling/data/sequence.fasta --outputDir /home/vidhi/project_4th_sem/Docker_var_calling/data/results

workflow {
    params.reads = file(params.reads ?: 'data/*.fastq')  
    params.reference = file(params.reference ?: 'data/*.fasta') 
    params.outputDir = params.outputDir ?: 'results'
   

    qc_result = quality_control(params.reads)
    read1 = file("/home/vidhi/project_4th_sem/check/data/SRR12345678_1.fastq")
    read2 = file("/home/vidhi/project_4th_sem/check/data/SRR12345678_2.fastq")

    trimmed_result = trim_adapters(read1, read2)
    adapter_trimming()
    alignment()
    counting()
    normalization()
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

process alignment {
    input:
    path trimmed_files from adapter_trimming.out

    output:
    path "results/aligned/*.bam"

    script:
    """
    mkdir -p ${params.outputDir}/aligned/
    echo "Directory ${params.outputDir}/aligned/ exists"
    
    # need to specify the paired reads for trimming
    trimmed1=\$(echo ${trimmed_result} | cut -d' ' -f1)
    trimmed2=\$(echo ${trimmed_result} | cut -d' ' -f2)
    hisat2 -x genome_index -U -1 \${trimmed1} -2 \${trimmed2} | samtools view -bS - > results/aligned/${trimmed_files.simpleName}.bam
    """
}

process counting {
    input:
    path bam_files from alignment.out

    output:
    path "results/counts/*.txt"

    script:
    """
    mkdir -p results/counts
    featureCounts -a annotation.gtf -o results/counts/${bam_files.simpleName}_counts.txt ${bam_files}
    """
}

process normalization {
    input:
    path count_files from counting.out

    output:
    path "results/normalization/*.txt"

    script:
    """
    mkdir -p results/normalization
    Rscript -e '
    library(edgeR)
    counts <- read.delim("${count_files}", row.names=1)
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge, method="TMM")
    write.table(cpm(dge), "results/normalization/CPM_${count_files.simpleName}.txt", sep="\\t")
    write.table(rpkm(dge), "results/normalization/FPKM_${count_files.simpleName}.txt", sep="\\t")
    '
    """
}




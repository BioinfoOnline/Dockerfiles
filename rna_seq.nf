params.data = "/home/vidhi/project_4th_sem/rna_seq/data"
params.output = "/home/vidhi/project_4th_sem/rna_seq/results"
params.adapter = "AGATCGGAAGAGC"

process quality_control {
    input:
    path fastq_files from file("${params.data}/*.fastq")

    output:
    path "results/qc/*.fastqc.zip"
    path "results/qc/*.fastqc.html"

    script:
    """
    mkdir -p results/qc
    fastqc ${fastq_files} -o results/qc
    """
}

process adapter_trimming {
    input:
    path fastq_files from file("${params.data}/*.fastq")

    output:
    path "results/trimmed/*.fastq"

    script:
    """
    mkdir -p results/trimmed
    trimmomatic SE -phred33 ${fastq_files} results/trimmed/$(basename ${fastq_files} .fastq)_trimmed.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50

    """
}

process alignment {
    input:
    path trimmed_files from adapter_trimming.out

    output:
    path "results/aligned/*.bam"

    script:
    """
    mkdir -p results/aligned
    hisat2 -x genome_index -U ${trimmed_files} | samtools view -bS - > results/aligned/${trimmed_files.simpleName}.bam
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

workflow {
    quality_control()
    adapter_trimming()
    alignment()
    counting()
    normalization()
}


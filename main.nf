#!/usr/bin/env nextflow

params.read1 = "s3://averafastq/everything_else/NA18238-b_S9_10k_1.fastq.gz"
params.read2 = "s3://averafastq/everything_else/NA18238-b_S9_10k_2.fastq.gz"
params.index = "s3://averagenomedb/Homo_sapiens/UCSC/hg19/star_genome"
params.gtf = "s3://averagenomedb/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
params.out = "s3://averatest/star_test/"
params.staropts = "--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical"

log.info "TRANSCRIPTOME QUANT P I P E L I N E"
log.info "================================="
log.info "Index              : ${params.index}"
log.info "Read 1             : ${params.read1}"
log.info "Read 2             : ${params.read2}"
log.info "Annotation         : ${params.gtf}"
log.info ""
log.info "Current home       : $HOME"
log.info "Current user       : $USER"
log.info "Current path       : $PWD"
log.info "Script dir         : $baseDir"
log.info "Working dir        : $workDir"
log.info "Output dir         : ${params.out}"
log.info ""

genome_index = file(params.index)
gft = file(params.gtf)
out = file(params.out)

process star {

    input:
    file genome_index
    file gtf
    file read1 from params.read1
    file read2 from params.read2
    file(out)
    
    output:
    file '*.Aligned.out.sorted.bam' into results
    file '*.Log.final.out' into results
    file '*.Log.out' into results
    file '*.Log.progress.out' into results
    file '*.SJ.out.tab' into results
 
    """
    prefix=\$(echo $read1 | sed 's/_.*//')
    cpu=\$(nproc)
    STAR --genome-dir $genome_index \\
    	 --sjdbGTFfile $gtf \\
    	 --readFilesIn $read1 $read2 \\
    	 --readFilesCommand zcat \\
    	 $params.staropts \\
    	 --runThreadN cpu \\ 
    	 --outFileNamePrefix $prefix
    """
}

results.subscribe { 
    log.info "Copying results to file: ${out}/${it.name}"
    it.copyTo(out)
 }
#!/usr/bin/env nextflow

params.read1 = Channel.fromPath("s3://averafastq/everything_else/*_1.fastq.gz")
params.read2 = Channel.fromPath("s3://averafastq/everything_else/*_2.fastq.gz")
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

gtf = file(params.gtf)
index = file(params.index)

process star {

    input:
    file gtf
    file read1 from params.read1
    file read2 from params.read2
    file index

    output:
    file '*.Aligned.out.sorted.bam' into results
    file '*.Log.final.out' into results
    file '*.Log.out' into results
    file '*.Log.progress.out' into results
    file '*.SJ.out.tab' into results

    """
    prefix=\$(echo $read1 | sed 's/_.*//')
    cpu=\$(grep -c "processor" /proc/cpuinfo)
    STAR --genomeDir $index \\
         --sjdbGTFfile $gtf \\
         --readFilesIn $read1 $read2 \\
         --readFilesCommand zcat \\
         $params.staropts \\
         --runThreadN \$cpu \\
         --outFileNamePrefix \$prefix
    """
}

results.subscribe {
    log.info "Copying results to file: ${params.out}/${it.name}"
    it.copyTo(out)
 }

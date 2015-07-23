#!/usr/bin/env nextflow

params.read1 = Channel.fromPath("s3://averafastq/everything_else/*_1.fastq.gz")
params.read2 = Channel.fromPath("s3://averafastq/everything_else/*_2.fastq.gz")
params.index = "/shared/Homo_sapiens/UCSC/hg19/star_2.4.2a_genome"
params.gtf = "s3://averagenomedb/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
params.out = "s3://averatest/star_test/"
params.towpassMode = "Basic"
params.quantMode = "TranscriptomeSAM GeneCount"
params.alignIntronMax = "200000"
params.alignMatesGapMax = "200000"
params.outFilterMismatchNoverLmax = "0.04"
params.outFilterMismatchNoverLmax = "RemoveNoncanonical"
params.outSAMtype = "BAM SortedByCoordinate"
params.outSAMunmapped = "Within"
params.outSAMattrRGline= "ID:$rgid SM:$rgsm PL:$rgpl LB:$rglb"
params.outSAMstrandField = "intronMotif"
params.chimSegmentMin = "25"
params.chimJunctionOverhangMin = "25"
params.outTmpDir = "/tmp/tmp"

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
    file '*.Aligned.sortedByCoord.out.bam' into bam
    file '*.Log.final.out' into results
    file '*.Log.out' into results
    file '*.Log.progress.out' into results
    file '*.SJ.out.tab' into results
    file '*.Chimeric.out.junction' into results
    file '*.Chimeric.out.sam' into results 
    file '*.ReadsPerGene.out.tab' into results
    file '*.Aligned.toTranscriptome.out.bam' into bam

    """
    prefix=\$(echo $read1 | sed 's/_.*//')
    cpu=\$(grep -c "processor" /proc/cpuinfo)
    STAR --genomeDir $index \\
         --sjdbGTFfile $gtf \\
         --readFilesIn $read1 $read2 \\
         --readFilesCommand cat \\
         --runThreadN \$cpu \\
         --outFileNamePrefix \$prefix
         --genomeDir $genomeDir \
		 --twopassMode $params.twopassMode \
		 --quantMode $params.quantMode \	
		 --alignIntronMax $params.alignIntronMax \
		 --alignMatesGapMax $params.alignMatesGapMax \
		 --outFilterMismatchNoverLmax $params.outFilterMismatchNoverLmax \
		 --outFilterIntronMotifs $params.outFilterIntronMotifs \
		 --outSAMtype $params.outSAMtype \
		 --outSAMunmapped $params.outSAMunmapped \
		 --outSAMattrRGline $params.outSAMattrRGline \
		 --outSAMstrandField $params.outSAMstrandField \
		 --chimSegmentMin $params.chimSegmentMin \
		 --chimJunctionOverhangMin $params.chimJunctionOverhangMin \
		 --outTmpDir $params.scratch
    """
}

process index {
	
	input:
	file '*.Aligned.sortedByCoord.out.bam' from bam
	file '*.Aligned.toTranscriptome.out.bam' from bam

	output:
	file '*.Aligned.sortedByCoord.out.bam.bai'
	file '*.Aligned.toTranscriptome.out.bam.bai'

	"""
	cpu=\$(grep -c "processor" /proc/cpuinfo)
	sambamba index -t  $\cpu -p
	"""
}

results.subscribe {
    log.info "Copying results to file: ${params.out}/${it.name}"
    it.copyTo(out)
}


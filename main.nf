params.run = 'runs.txt'
params.outdir = 'results'

process DOWNLOAD_SRA {
	tag "$sample"
	label 'download'

	input:
	val sample

	output:
	tuple val(sample), path("${sample}")

	script:

	"""

	prefetch ${sample} --max-size 50G --output-directory .

	"""

}

process SRA_TO_FASTQ {

	tag "$sample"
	label 'download'	

	publishDir "${params.outdir}/fastq", mode: 'copy'

	input:
	tuple val(sample), path(sra_dir)

	output:
	tuple val(sample),
		path("${sample}_1.fastq.gz"),
		path("${sample}_2.fastq.gz")

	script:
	"""
	mkdir -p tmp
	
	fasterq-dump ${sra_dir} \\
		--split-files \\
		--threads ${task.cpus} \\
		--temp tmp \\
		--outdir .

	gzip ${sample}_1.fastq ${sample}_2.fastq

	"""

}	

process FASTQC {

	tag "$sample"

	publishDir "${params.outdir}/fastqc", mode: 'copy'

	input:
	tuple val(sample), path(r1), path(r2)

	output:
	path("*_fastqc.html")
	path("*_fastqc.zip")

	script:

	"""

	fastqc ${r1} ${r2} \\
	--threads ${task.cpus} \\
	--outdir .

	"""
}

process INDEX_REFERENCE{

	input:
	path fasta
	output:
	path 'refidx"
	script:

	"""
	mkdir refidx
	cp '${fasta}' refidx/reference.fa
	bwa index refidx/reference.fa
	samtools faidx refidx/reference.fa
	"""
}


process ALIGN{

	tag "$sample"
	label 'align'
	publishDir "${params.outdir}/bam", mode:'copy'
	input:
	tuple val(sample), path(r1), path(r2)
	path refidx
	output:
	tuple val(sample),
		path("${sample}.bam"),
		path("${sample}.bam.bai"),
		emit:bams

	path "${sample}.flagstat.txt", emit: stats
	path "${sample}.markdup.txt", emit: duplicates
	script:

	"""
	bwa mem -t ${task.cpus} \\
		-R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' \\
		refidx/reference.fa '${r1}' '${r2}' \\
		| samtools sort -n -m 256M -o names.bam -
	samtools fixmate -m names.bam mates.bam
	samtools sort -m 256M -o sorted bam mates.bam

	samtools markdup -s \\
		-f '${sample}.markdup.txt \\
		sorted.bam '${sample}.bam' \\

	samtools index '${sample}.bam'

	samtools flagstat '${sample}.bam' \\ 
		> '${sample}.flagstat.txt	
	
	"""
}

process CALL_VARIANTS{

	label 'call'

	publishDir "$params.outdir/vcf", mode: 'copy'

	input:
	tuple path(bam_files), path(bai_files)
	path refidx

	output:
	tuple path('cohort.vcf.gz'),
		path('cohort.vcf.gz.csi'),
		emit: vcf

	path 'cohort.stats.txt', emit: stats

	script:
	def bams = bam_files.collect {"'${it}'"}.join(' ')

	"""
    	bcftools mpileup \\
        	--fasta-ref refidx/reference.fa \\
        	--min-MQ ${params.min_mapq} \\
        	--min-BQ ${params.min_baseq} \\
        	--max-depth ${params.max_depth} \\
        	--annotate FORMAT/DP,FORMAT/AD \\
        	-Ou ${bams} \\
        	| bcftools call \\
            		--ploidy 1 \\
            		--multiallelic-caller \\
            		--variants-only \\
            		-Ou \\
        	| bcftools norm \\
            		--fasta-ref refidx/reference.fa \\
            		-Ou \\
        	| bcftools sort \\
            		--max-mem 256M \\
            		-Oz \\
            		-o cohort.vcf.gz

    bcftools index cohort.vcf.gz

    bcftools stats cohort.vcf.gz > cohort.stats.txt
    """

}

process BUILD_SNPEFF{

	label 'annotation'

	input:
	path fasta
	path gff

	output:
	path 'annotation_db'

	script:

	def heapMb = (task.memory.toMega() *0.75) as int

	"""

	mkdir -p annotation_db/data/project

	cp '${fasta}' annotation_db/data/project/sequences.fa
	cp '${gff}' annotation_db/data/project/genes.gff

	printf 'data.dir = ./data\\nproject.genome : Project_reference\\n' \\
        > annotation_db/snpEff.config

    	snpEff -Xms128m -Xmx${heapMb}m build \\
        	-gff3 \\
        	-noCheckCds \\
        	-noCheckProtein \\
        	-c annotation_db/snpEff.config \\
        	project
   	 """

}

process ANNOTATE{

	label 'annotation'
	
	publishDir "${params.outdir}/annotation", mode: 'copy'

	input:

	tuple path(vcf), path(csi)
	path annotation_db

	output:

	tuple path('cohort.annotated.vcf.gz'),
		path('cohort.annotated.vcf.gz.csi')
		emit: vcf

	path 'snpeff_summary.html', emit: report
		optional: true,
		emit: genes

	script:
	def heapMb = (task.memory.toMega()*0.75) as int

	"""

	snpEff -Xms128m -Xmx${heapMb}m ann \\
        	-noDownload \\
        	-c annotation_db/snpEff.config \\
       		-stats snpeff_summary.html \\
       		project '${vcf}' \\
        	| bcftools view -Oz -o cohort.annotated.vcf.gz

    bcftools index cohort.annotated.vcf.gz	

	"""

}

process UNIQUE_VARIANTS{}


workflow {
	
	runs = Channel
	.fromPath(params.runs, checkIfExists: true)
	.splitText()
	.map {it.trim()}
	.filter{ it && !it.startsWith('#') }

	DOWNLOAD_SRA(runs)
	SRA_TO_FASTQ(DOWNLOAD_SRA.out)
	FASTQC(SRA_TO_FASTQ.out)

}


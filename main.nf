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
	
	fasterq-dump ${sra_dir} \
		--split-files \
		--threads ${task.cpus} \
		--temp tmp \
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

	fastqc ${r1} ${r2} \
	--threads ${task.cpus} \
	--outdir .

	"""
}

process INDEX_REFERENCE{}
process ALIGN{}
process CALL_VARIANTS{}
process BUILD_SNPEFF{}
process ANNOTATE{}
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


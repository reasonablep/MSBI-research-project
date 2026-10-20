params.run = 'runs.txt'
params.outdir = 'results'

process DOWNLOAD_SRA {
	tag "$run"

	cpus 1
	memory '4 GB'
	time '4h'

	input:
	val run

	output:
	tuple val(run), path("${run}")

	script:

	"""

	prefetch ${run} --max-size 50G --output-directory .

	"""

}

process SRA_TO_FASTQ {

	tag "$run"
	

	cpus 1
	memory '8GB'
	time '3h'

	publishDir "${params.outdir}/fastq", mode: 'copy'

	input:
	tuple val(run), path(sra_dir)

	output:
	tuple val(run),
		path("${run}_1.fastq.gz"),
		path("${run}_2.fastq.gz")

	script:
	"""
	mkdir -p tmp
	
	fasterq-dump ${sra_dir} \
		--split-files \
		--threads ${task.cpus} \
		--temp tmp \
		--outdir .

	gzip ${run}_1.fastq ${run}_2.fastq

	"""

}	

process FASTQC {

	tag "$run"

	cpus 1
	memory '4 GB'
	time '1h'

	publishDir "${params.outdir}/fastqc", mode: 'copy'

	input:
	tuple val(run), path(read1), path(read2)

	output:
	path("*_fastqc.html")
	path("*_fastqc.zip")

	script:

	"""

	fastqc ${read1} ${read2} \
	--threads ${task.cpus} \
	--outdir .

	"""
}


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



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
	cp "${fasta}" refidx/reference.fa
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
		-R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' refidx/reference.fa '${r1}' '${r2}' samtools sort -n -m 256M -o names.bam -
	samtools fixmate -m names.bam mates.bam
	samtools sort -m 256M -o sorted.bam mates.bam

	samtools markdup -s -f '${sample}.markdup.txt sorted.bam '${sample}.bam'

	samtools index '${sample}.bam'

	samtools flagstat '${sample}.bam' > '${sample}.flagstat.txt'	
	
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

process CALL_FREEBAYES{

	label 'call'

	publishDir "{params.outdir}/freebayes" mode: 'copy'
	
	input:
	tuple path(bam_files), path(bai_files)
	path refidx

	output:
	tuple path('freebayes.vcf.gz'),
		path('freebayes.vcf.gz.csi')
		emit: vcf

	 
	script:
	def bams= bam_files.collect{ "'${it}'"  }.join(' ')

	"""
	
	freebayes -f refidx/reference.fa -p 1 ${bams} \
		| bcftools norm -f refidx/reference.fa -Ou \
		| bcftools sort -Oz -o freebayes.vcf.gz
		| bcftools index --csi freebayes.vcf.gz

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

process UNIQUE_VARIANTS{

	publishDir "${params.outdir}/candidates", mode: 'copy'

	input:
	tuple path(vcf), path(csi)

	output:
	path '*.tsv'

	script:

	"""
	
	unique_variants.py \\
		--vcf '{vcf}'
		--min-depth ${params.min_depth} \\
		--min-qual ${params.min_qual} \\
		--min-alt-fraction ${params.min_alt_fraction} \\
		--max-ref-alt-fraction ${params.max_ref_alt_fraction}
	"""

}


workflow {
	if (params.help){

		log.info 'See README.md. Local test: nextflow run main.nf -profile local,test'
	} else {
		if(!params.reference){
			error 'Provide --reference references/genome.fa(uncompressed FASTA)'}
		if ((params.reads ? 1 : 0) + (params.runs ? 1 : 0) != 1) {
			error 'Provide exactly one --reads or --runs'
		}
		
		if(
			params.min_depth <1 ||
			params.min_qual <0 ||
			params.min_mapq < 0 ||
			params.min_baseq <0 
	){
			error 'Depth must be positive; quality thresholds non-negative'
	}
	
		if(
			params.min_alt_fraction <=0.5 ||
			params.min_alt_fraction > 1 ||
			params.max_ref_alt_fraction <0 ||
			params.max_ref_alt_fraction >= params.min_alt_fraction
	){
			error 'Invalid allele fraction thresholds"'
	}
		fasta = Channel.value(
		file(params.reference, checkIfExists: true)
		)

		if (params.reads){

			reads_ch = Channel
				.fromFilePairs(
					params.reads,
					size: -1,
					checkIfExists:true
				)	
				.map{sample, pair ->
					if(!(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)){
						error "Unsafe sample ID: ${sample}"
					}
					
					if(pair.size()!=2){
						error "Expected exactly two mates for ${sample}; found ${pair.size()}"	
					}

					tuple(sample,pair[0],pair[1])

				}
				
			} else{

				runs_ch=Channel
					.fromPath(params.runs, checkIfExists: true)
					.splitText()
					.map{it.trim()}
					.filter{it && !it.startsWith('#')}
					.map{ run-> 
						if(!(run==~ /[SED]RR[0-9]+/)){
							error "Invalid accession: ${run}"
						}

						run
					
					}

					.unique()
					.ifEmpty{
						error 'No accessions in --runs file'
					}
				DOWNLOAD_SRA(runs_ch)
				SRA_TO_FASTQ(DOWNLOAD_SRA.out)

				reads_ch = SRA_TO_FASTQ.out
						
			}

			FASTQC(read_ch)
			INDEX_REFERENCE(fasta)
			ALIGN(reads_ch, INDEX_REFERENCE.out)
			
			cohort = ALIGN.out.bams
				.toList()
				.map{ rows->
					
					if(rows.size() <2){
					error 'At least two samples required for strain-unique comparison
					}

					def ordered = rows.sort {a, b -> 
						a[0] <=> b[0]
					}

					tuple(
						ordered.collect { it[1] },
						ordered.collect { it[2] }
					)

				}
			CALL_VARIANTS(cohort, INDEX_REFERENCE.out)
			CALL_FREEBAYES(cohort, INDEX_REFERENCE.out)			
			if (params.gff){
				BUILD_SNPEFF(
					fasta,
					Channel.value(
						file(params.gff, checkIfExists: true)
					)
				)

				ANNOTATE(
					CALL_VARIANTS.out.vcf,
					BUILD_SNPEFF.out
				)

				UNIQUE_VARIANTS(ANNOTATE.out.vcf)

			} else {
				UNIQUE_VARIANTS(CALL_VARIANTS.out.vcf)
			}
				 
		}
	}		


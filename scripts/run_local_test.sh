#!/bin/bash

set -euo pipefail

# Toggle annotations on or off

cd "$(dirname"$0")/.."
annotate = 0
if [ "${1:-}" = '--annotation' ]; then
	command -v snpEff >/dev/null 2>&1 || { echo 'Missing snpEff' >&2; exit 1; }
	annotate=1

elif [ "$#" -gt 0 ]; then
	echo 'Usage: bash scripts/run_local_test.sh [--annotation]' >&2
	exit 1

fi

for tool in python3 java nextflow fastqc bwa samtools bcftools; do

	if ! command -v "$tool" >/dev/null 2>&1; then
	printf 'Missing tool: %s. See environment.yml and README.md.\n' "$tool" >&2
	exit 1	
	fi
done
if [ ! -d testdata ]; then
	python3 scripts/make_test_data.py
fi

# Make report folder after path checks verified 

mkdir -p reports/local

# Build nextflow command, execute, check results

set -- nextflow run main.nf -profile local,test -resume
if [ "$annotate" = '1' ]; then
    set -- "$@" --gff testdata/references/toy.gff3
fi
if [ "${NSPORA_REPORTS:-1}" = '1' ]; then
    set -- "$@" -with-report reports/local/report.html \
        -with-trace reports/local/trace.tsv -with-timeline reports/local/timeline.html
fi
"$@"
if [ "$annotate" = '1' ]; then
    python3 scripts/check_test_results.py --require-annotation
else
    python3 scripts/check_test_results.py
fi

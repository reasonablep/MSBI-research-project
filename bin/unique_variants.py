#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
from pathlib import Path

# Examine one sample at one variant position adn classify evidence
def inspect_sample(fmt, text, args):
	fields = dict(zip(fmt, text.split(':')))
	try:
		gt = fields['GT']
		dp = int(fields['DP'])
		ad = [int(v) for v in fields['AD'].split(',')]
		if gt not in ('0', '1') or len(ad) !=2 or min(ad) <0:
			return 'unknown', dp, None
		total = sum(ad)
		if dp < args.min_depth or total < args.min_depth:
			return 'unknown', dp, None
		fraction = ad[1]/total
		if gt =='1' and fraction >= args.min_alt_fraction:
			return 'alt', dp, fraction
		if gt=='0' and fraction <= args.max_ref_alt_fraction:
			return 'ref', dp, fraction
		return 'mixed', dp, fraction
	except (KeyError, ValueError, ZeroDivisionError):
		return 'unknown', None, None

# Extract annotations of positions that qualify as unique variants
def annotations(info, alt):
    raw= next((v[4:] for v in info.split(';') if v.startswith('ANN=')), '')
    entries = [v.split('|') for v in raw.split(',') if v]
    entries= [v for v in entries if len(v) > 4 and v[0]== alt]
    return [','.join(sorted({v[i] for v in entries if v[i]})) for i in (3, 1, 2)]+ [raw]

# Coordinate the comparison by reading the VCF, calling inspect_sample, decide uniqueness, write results

def main():
    p= argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf', required= True)
    p.add_argument('--outdir', default='.')
    p.add_argument('--min-depth', type= int, default=10)
    p.add_argument('--min-qual', type= float, default=30)
    p.add_argument('--min-alt-fraction', type=float, default=0.8)
    p.add_argument('--max-ref-alt-fraction', type=float, default=0.1)
    args = p.parse_args()

    if args.min_depth < 1 or args.min_qual < 0 or not (
        0 <= args.max_ref_alt_fraction <args.min_alt_fraction <= 1
        and args.min_alt_fraction >0.5
    ):
        p.error('Invalid quality/depth/allele-fraction thresholds')

#Prepare the output directory and input reader

    out =Path(args.outdir)
    out.mkdir(parents= True, exist_ok=True)
    opener= gzip.open if args.vcf.endswith('.gz') else open
    header= ['sample', 'chrom', 'pos', 'ref', 'alt', 'qual', 'depth',
              'alt_fraction', 'gene', 'consequence', 'impact', 'ANN']
    handles= []
    samples= None

# Open the input VCF and combined output files in tsv formt

    try:
        with opener(args.vcf, 'rt') as vcf, (out / 'unique_variants.tsv').open('w') as unique_f, \
                (out / 'variant_audit.tsv').open('w') as audit_f:
            unique= csv.writer(unique_f, delimiter='\t')
            unique.writerow(header)
            audit= csv.writer(audit_f, delimiter='\t')
            audit.writerow(['chrom', 'pos', 'ref', 'alt', 'decision', 'sample_states'])
            writers={}

# Read sample names from the VCF header

            for line in vcf:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    samples= line.rstrip('\n').split('\t')[9:]
                    if len(samples)< 2 or len(set(samples))!= len(samples):
                        raise ValueError('Need at least two distinct jointly genotyped samples')
                    for sample in samples:

                        if not sample or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-' for c in sample):
                            raise ValueError('Invalid sample ID in VCF')
                        f= (out / (sample + '.unique.tsv')).open('w')
                        handles.append(f)
                        writers[sample]= csv.writer(f, delimiter='\t')
                        writers[sample].writerow(header)
                    continue
                if line.startswith('#') or not line.strip():
                    continue
                if samples is None:
                    raise ValueError('Missing VCF sample header')

# Separate variant record into fields

                row = line.rstrip('\n').split('\t')
                if len(row) != 9 + len(samples):
                    raise ValueError('Wrong number of sample columns')
                chrom, pos, _, ref, alt, qual, filt, info, format_string = row[:9]

# Inspect the samples and decides if the variant is shared or unique
                states= {}
                if ',' in alt or alt == '*' or '<' in alt:
                    decision= 'complex_or_multiallelic_review'
                elif filt not in ('.', 'PASS') or qual == '.' or float(qual) < args.min_qual:
                    decision= 'site_quality_fail'
                else:
                    states ={s: inspect_sample(format_string.split(':'), value, args)
                              for s, value in zip(samples, row[9:])}
                    labels =[value[0] for value in states.values()]
                    if 'unknown' in labels:
                        decision= 'unassessable_in_one_or_more_samples'
                    elif 'mixed' in labels:
                        decision= 'mixed_allele_evidence'
                    elif labels.count('alt') > 1:
                        decision= 'shared_variant'
                    elif labels.count('alt') == 1 and labels.count('ref') == len(samples) - 1:
                        decision= 'unique_within_cohort'

# Record accepted unique variants and finds the source sample, recording sample depth

                        sample= next(s for s, value in states.items() if value[0] == 'alt')
                        _, dp, fraction = states[sample]
                        result = [sample, chrom, pos, ref, alt, qual, dp, f'{fraction:.4f}']
                        result+= annotations(info, alt)
                        unique.writerow(result)
                        writers[sample].writerow(result)
                    else:
                        decision = 'no_confident_alternate'

# Record every decision 

                audit.writerow([chrom, pos, ref, alt, decision,
                                json.dumps(states, sort_keys=True, separators=(',', ':'))])
            if samples is None:
                raise ValueError('Missing VCF sample header')

# Closes per-sample files and start the program via outside execution

    finally:
        for handle in handles:
            handle.close()

if __name__=='__main__':
	main()


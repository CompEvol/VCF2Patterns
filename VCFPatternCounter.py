# Author: Walter Xie
# Last update: 5 Apr 2016

import vcf
import sys
import re
from Bio import SeqIO
from collections import Counter

# return a tuple of sample names on 2 copies
def get_samples(vcf_reader):
	samples = ()
	#print (vcf_reader.samples)
	for sn in vcf_reader.samples:
		s1 = sn + "_1"
		s2 = sn + "_2"
		samples += (s1,s2,)
	return samples

# return a string of one-site pattern in a site of all samples on 2 copies
# if ALT contains insertion/deletion then return empty string
def get_one_site_patttern(record, ref, alt, verbose=False):
	site_patt=""
	if len(ref) != 1:
		return ""

	for sample in record.samples:
		#sample = record.samples[1]
		gt = sample['GT'].strip()

		if verbose:
			print(gt)

		if gt == ".":
			site_patt += ref * 2
		else:
			# ignore phase/unphase
			for gt2 in re.split('\||/', gt):
				gt2 = int(gt2)
				if gt2 > 0:
					alt_gt = str(alt[gt2-1])
				elif gt2 == 0:
					alt_gt = ref
				else:
					raise Exception("Invaild genotype, sample = {}, GT = {} !"
						.format(sample.sample, gt))
				site_patt += alt_gt

			if len(alt_gt) > 1:
				return ""

	if verbose:
		print(record)

	if len(site_patt) != len(record.samples)*2:
		print(site_patt)
		raise Exception("One site pattern has incorrect number in {} "
			"position {} of sample {}, where {} != {} !"
			.format(record.CHROM, record.POS, sample.sample, len(site_patt), len(record.samples)*2))

	return site_patt

# count all bases in reference sequences
# return a dictionary of nt => counts
def count_ref_bases(reference_sequences):
	ref_bases = Counter()
	for ref_seq in reference_sequences:
		#ref_seq = next(reference_sequences)
		seq = ref_seq.seq
		tot_nt = 0
		for nt in ['A', 'C', 'G', 'T', 'N']:
			c_nt = seq.count(nt)
			tot_nt += c_nt
			ref_bases[nt] += c_nt
		if len(seq) != tot_nt:
			raise Exception("There are {} nucleotide bases, but find {} in total !"
				.format(len(seq), tot_nt))

	return ref_bases

######### main

#vcf_f_n = sys.argv[1]
vcf_f_n="./adelie/allmodern.24.allancient.22.q20.vcf"
#ref_f_n = sys.argv[2]
ref_f_n="./adelie/adelie.allscaffolds.fa"
msg_freq = 50000

vcf_reader = vcf.Reader(open(vcf_f_n, 'r'))
print("Count patterns from", vcf_f_n)

rows_not_count_file = open("RowsNotCount.txt", "w")
patterns_file = open("Patterns.txt", "w")

samples=get_samples(vcf_reader)

# bases are taken for non-constant sites, nt in ref => counts
ref_bases_taken = Counter()
# pattern => counts
patterns = Counter()
row = 0
for record in vcf_reader:
	#record = next(vcf_reader)
	#print (record)

	ref = record.REF
	alt = record.ALT
	#chrom = record.CHROM
	#pos = record.POS

	site_patt = get_one_site_patttern(record, ref, alt)
	# only pick up mutation, ignore insertion/deletion
	if len(site_patt) > 0:
		ref = ref.upper()
		ref_bases_taken[ref] += 1
		patterns[site_patt] += 1
	else:
		print(record, file=rows_not_count_file)
	row += 1

	if row % msg_freq == 0:
		print("Find", len(patterns), "patterns from", row, "rows.")

print("Find", len(patterns), "patterns from total", row, "rows.")

# print patterns
print('\t'.join(samples), file=patterns_file)
for k, v in patterns.items():
	print('{}\t{}'.format(k, v), file=patterns_file)

# count all constant sites and add to patterns
reference_sequences = SeqIO.parse(open(ref_f_n),'fasta')
print("There are", len(reference_sequences), "sequences in the reference", ref_f_n)

ref_bases = count_ref_bases(reference_sequences)

diff = set(ref_bases.keys()) - set(ref_bases_taken.keys())
if len(diff) > 0:
	print("Warning: find different nucleotide bases between reference and VCF, "
			"diff = {} !".format(diff))

for k in ref_bases.keys():
	constant_site = k * len(samples)
	if k not in ref_bases_taken:
		print("Warning: cannot find nucleotide base {}, skip it !".format(k))
		continue
	count = ref_bases[k] - ref_bases_taken[k]
	if count < 1:
		print("Warning: constant site having all {}s does not have a positive count, "
				"count = {} !".format(k, count))
		continue
	print('{}\t{}'.format(constant_site, count), file=patterns_file)

patterns_file.close()
rows_not_count_file.close()



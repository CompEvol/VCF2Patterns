# Author: Walter Xie
# Last update: 5 Apr 2016

import vcf
import sys
import re
from Bio import SeqIO
from collections import Counter
import json

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
		return "", 0

	n_samples = 0
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
					n_samples += 1
				elif gt2 == 0:
					alt_gt = ref
				else:
					raise Exception("Invaild genotype, sample = {}, GT = {} !"
						.format(sample.sample, gt))
				site_patt += alt_gt

			if len(alt_gt) > 1:
				return "", 0
			
	if verbose:
		print(record)

	if len(site_patt) != len(record.samples)*2:
		print(site_patt)
		raise Exception("One site pattern has incorrect number in {} "
			"position {} of sample {}, where {} != {} !"
			.format(record.CHROM, record.POS, sample.sample, len(site_patt), len(record.samples)*2))
	
	if n_samples > len(record.samples)*2 and n_samples < 1:
		raise Exception("Incorrect number of indels {} should <= total number "
			"of samples {}, {} position {} of sample {} !"
			.format(n_samples, len(record.samples)*2), record.CHROM, record.POS, sample.sample)
	
	return site_patt, n_samples

def get_indel_count(record, ref, alt, verbose=False):
	n_samples = 0
	
	for sample in record.samples:
		#sample = record.samples[1]
		gt = sample['GT'].strip()

		if verbose:
			print(gt)

		if gt != ".":
			# ignore phase/unphase
			for gt2 in re.split('\||/', gt):
				gt2 = int(gt2)
				if gt2 > 0:
					n_samples += 1
					
	if n_samples > len(record.samples)*2 and n_samples < 1:
		raise Exception("Incorrect number of indels {} should <= total number "
			"of samples {}, {} position {} of sample {} !"
			.format(n_samples, len(record.samples)*2), record.CHROM, record.POS, sample.sample)
	
	return n_samples


# count all bases in reference sequences
# return a dictionary of nt => counts
def count_ref_bases(reference_sequences):
	n_seq = 0
	for ref_seq in reference_sequences:
		#ref_seq = next(reference_sequences)
		seq = ref_seq.seq
		ref_bases = count_bases(seq)
		n_seq += 1
						
	print("There are {} reference sequences in total".format(n_seq))	
	return ref_bases

def count_bases(seq):
	ref_bases = Counter()
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
vcf_f_n="./adelie/allmodern.24.8x.2.allancient.22.q20.vcf"
#ref_f_n = sys.argv[2]
ref_f_n="./adelie/adelie.allscaffolds.fa"
msg_freq = 100000

vcf_reader = vcf.Reader(open(vcf_f_n, 'r'))
print("Count patterns from", vcf_f_n)

rows_not_count_file = open("RowsNotCount.txt", "w")

samples=get_samples(vcf_reader)

# bases are taken for non-constant sites including ins/del, nt in ref => counts
ref_bases_taken = Counter()
# SNP patterns, pattern => counts
patterns = Counter()
# number of samples mutated on the site, number of samples => counts
mu_frequencies = Counter()
indel_frequencies = Counter()
row = 0
n_c_row = 0
for record in vcf_reader:
	#record = next(vcf_reader)
	#print (record)

	ref = record.REF
	alt = record.ALT
	#chrom = record.CHROM
	#pos = record.POS

	# all changes should be taken out from reference including ins/del
	ref_bases_taken += count_bases(ref)
	
	site_patt, mu_freq = get_one_site_patttern(record, ref, alt)
	# only pick up mutation, ignore insertion/deletion
	if len(site_patt) > 0:
		ref = ref.upper()
		patterns[site_patt] += 1
		mu_frequencies[mu_freq] +=1
	else:
		n_c_row += 1
		n_c_freq = get_indel_count(record, ref, alt)
		indel_frequencies[n_c_freq] += 1
		print(record, file=rows_not_count_file)
		
	row += 1

	if row % msg_freq == 0:
		print("Find {} patterns from {} rows.".format(len(patterns), row))

print("Find {} patterns from total {} rows, where {} rows are ignored.".format(len(patterns), row, n_c_row))

with open('patterns.json', 'w') as fp:
    json.dump(patterns, fp)
with open('ref_bases_taken.json', 'w') as fp:
    json.dump(ref_bases_taken, fp)

with open('indel_frequencies.txt', 'w') as fp:
	for k, v in indel_frequencies.items():
		print('{}\t{}'.format(k, v), file=fp)

indel_freq_sum = sum(indel_frequencies.values())
if n_c_row != indel_freq_sum:
	raise Exception("Incorrect indel frequencies {}, it should = {} !".format(indel_freq_sum, n_c_row))

with open('mu_frequencies.txt', 'w') as fp:
	for k, v in mu_frequencies.items():
		print('{}\t{}'.format(k, v), file=fp)

mu_freq_sum = sum(mu_frequencies.values())
if (row - n_c_row) != mu_freq_sum:
	raise Exception("Incorrect mutation frequencies {}, it should = {} !".format(mu_freq_sum, (row - n_c_row)))

# print patterns
patterns_file = open("Patterns.txt", "w")

print('\t'.join(samples), file=patterns_file)
for k, v in patterns.items():
	print('{}\t{}'.format(k, v), file=patterns_file)

# count all constant sites and add to patterns
reference_sequences = SeqIO.parse(open(ref_f_n),'fasta')
print("Loading reference sequences from", ref_f_n)

ref_bases = count_ref_bases(reference_sequences)
with open('ref_bases.json', 'w') as fp:
    json.dump(ref_bases, fp)

diff = set(ref_bases.keys()) - set(ref_bases_taken.keys())
if len(diff) > 0:
	print("Warning: find different nucleotide bases between reference and VCF ref, "
			"diff = {} !".format(diff))

for k in ref_bases.keys():
	constant_site = k * len(samples)
	if k not in ref_bases_taken:
		print("Warning: cannot find nucleotide base {} in dict ref_bases_taken, skip it !".format(k))
		continue
	count = ref_bases[k] - ref_bases_taken[k]
	if count < 1:
		print("Warning: constant site having all {}s does not have a positive count, "
				"count = {} !".format(k, count))
		continue
	print('{}\t{}'.format(constant_site, count), file=patterns_file)

patterns_file.close()
rows_not_count_file.close()

# count total constant site bases
constant_site_bases = sum(ref_bases.values()) - sum(ref_bases_taken.values())
print("Total constant site bases =", constant_site_bases) 




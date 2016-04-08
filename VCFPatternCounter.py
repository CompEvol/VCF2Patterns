# Author: Walter Xie
# Last update: 5 Apr 2016

import vcf
import sys
import re
from Bio import SeqIO
from collections import Counter
import json

# return a tuple of sample names on 2 copies
def get_samples(vcf_reader, ind_age):
	samples = ()
	#print (vcf_reader.samples)
	for sn in vcf_reader.samples:
		if ind_age == 0 or (ind_age == 1 and len(sn) <= ancient_s_n_len) or (ind_age == 2 and len(sn) > ancient_s_n_len):
			s1 = sn + "_1"
			s2 = sn + "_2"
			samples += (s1,s2,)
	return samples

# return a string of one-site pattern in a site of all samples on 2 copies
# if ALT contains insertion/deletion then return empty string
def get_one_site_patttern(record, ref, alt, ind_age, verbose=False):
	site_patt=""
	if len(ref) != 1:
		return "", 0

	freq = 0
	n_sample = 0
	for sample in record.samples:
		#sample = record.samples[1]
		gt = sample['GT'].strip()
		sample_name = sample.sample
		
		if verbose:
			print(gt)

		if ind_age == 0 or (ind_age == 1 and len(sample_name) <= ancient_s_n_len) or (ind_age == 2 and len(sample_name) > ancient_s_n_len):
			n_sample += 1
			if gt == ".":
				site_patt += ref * 2
			else:
				# ignore phase/unphase
				for gt2 in re.split('\||/', gt):
					gt2 = int(gt2)
					if gt2 > 0:
						alt_gt = str(alt[gt2-1])
						freq += 1
					elif gt2 == 0:
						alt_gt = ref
					else:
						raise Exception("Invaild genotype, sample = {}, GT = {} !"
							.format(sample_name, gt))
					site_patt += alt_gt

				if len(alt_gt) > 1:
					return "", 0
			
	if verbose:
		print(record)

	if len(site_patt) != n_sample*2:
		print(site_patt)
		raise Exception("One site pattern has incorrect number in {} "
			"position {} of sample {}, where {} != {} !"
			.format(record.CHROM, record.POS, sample_name, len(site_patt), n_sample*2))
	
	if freq > n_sample*2 and freq < 1:
		raise Exception("Incorrect number of SNPs {} should <= total number "
			"of samples {}, {} position {} of sample {} !"
			.format(freq, n_sample*2, record.CHROM, record.POS, sample_name))
	
	return site_patt, freq

def get_indel_count(record, ref, alt, ind_age, verbose=False):
	freq = 0
	n_sample = 0
	for sample in record.samples:
		#sample = record.samples[1]
		gt = sample['GT'].strip()
		sample_name = sample.sample

		if verbose:
			print(gt)

		if ind_age == 0 or (ind_age == 1 and len(sample_name) <= ancient_s_n_len) or (ind_age == 2 and len(sample_name) > ancient_s_n_len):
			n_sample += 1
			if gt != ".":
				# ignore phase/unphase
				for gt2 in re.split('\||/', gt):
					gt2 = int(gt2)
					if gt2 > 0:
						freq += 1
					
	if freq > n_sample*2 and freq < 1:
		raise Exception("Incorrect number of indels {} should <= total number "
			"of samples {}, {} position {} of sample {} !"
			.format(freq, n_sample*2, record.CHROM, record.POS, sample_name))
	
	return freq


# count all bases in reference sequences
# return a dictionary of nt => counts
def count_ref_bases(reference_sequences):
	ref_bases = Counter()
	n_seq = 0
	nt_tot = 0
	for ref_seq in reference_sequences:
		#ref_seq = next(reference_sequences)
		seq = ref_seq.seq
		ref_bases += count_bases(seq)
		nt_tot += len(seq)
		n_seq += 1
	print("There are {} reference sequences and {} nucleotide bases in total."
		.format(n_seq, nt_tot))	
	
	ref_bases_sum = sum(ref_bases.values())
	if ref_bases_sum != nt_tot:
		raise Exception("Find {} nucleotide bases, but there are {} in reference sequences file !"
			.format(ref_bases_sum, nt_tot))

	return ref_bases

def count_bases(seq):
	bases = Counter()
	nt_tot = 0
	for nt in ['A', 'C', 'G', 'T', 'N']:
		c_nt = seq.count(nt)
		bases[nt] += c_nt
		nt_tot += c_nt
	if len(seq) != nt_tot:
		raise Exception("There are {} nucleotide bases, but find {} in total !"
			.format(len(seq), nt_tot))
	return bases

######### main

#vcf_f_n = sys.argv[1]
#vcf_f_n="./adelie/1klines.vcf"
vcf_f_n="./adelie/allmodern.24.8x.2.allancient.22.q20.vcf"
#ref_f_n = sys.argv[2]
ref_f_n="./adelie/adelie.allscaffolds.fa"
#ind_age = sys.argv[3]
ind_age = 0 # 0 chooses all, 1 chooses just the modern individuals, 2 chooses just ancient 

msg_freq = 100000
# ancient sample names len() > 9
ancient_s_n_len = 9

vcf_reader = vcf.Reader(open(vcf_f_n, 'r'))
print("Count patterns from", vcf_f_n)

rows_not_count_file = open("RowsNotCount.txt", "w")

samples=get_samples(vcf_reader, ind_age)

# bases are taken for non-constant sites including ins/del, nt in ref => counts
ref_bases_taken = Counter()
# SNP patterns, pattern => counts
patterns = Counter()
# number of samples mutated on the site, number of samples => counts
snp_frequencies = Counter()
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
	
	site_patt, snp_freq = get_one_site_patttern(record, ref, alt, ind_age)
	# only pick up mutation, ignore insertion/deletion
	if len(site_patt) > 0:
		ref = ref.upper()
		patterns[site_patt] += 1
		snp_frequencies[snp_freq] +=1
	else:
		n_c_row += 1
		n_c_freq = get_indel_count(record, ref, alt, ind_age)
		indel_frequencies[n_c_freq] += 1
		print(record, file=rows_not_count_file)
		
	row += 1

	if row % msg_freq == 0:
		print("Find {} patterns from {} rows.".format(len(patterns), row))

rows_not_count_file.close()
print("Find {} patterns from total {} rows, where {} rows are ignored.".format(len(patterns), row, n_c_row))

#with open('patterns.json', 'w') as fp:
#    json.dump(patterns, fp)
with open('ref_bases_taken.json', 'w') as fp:
    json.dump(ref_bases_taken, fp)

with open('indel_frequencies.txt', 'w') as fp:
	for k, v in indel_frequencies.items():
		print('{}\t{}'.format(k, v), file=fp)

indel_freq_sum = sum(indel_frequencies.values())
if n_c_row != indel_freq_sum:
	raise Exception("Incorrect indel frequencies {}, it should = {} !".format(indel_freq_sum, n_c_row))

with open('snp_frequencies.txt', 'w') as fp:
	for k, v in snp_frequencies.items():
		print('{}\t{}'.format(k, v), file=fp)

snp_freq_sum = sum(snp_frequencies.values())
if (row - n_c_row) != snp_freq_sum:
	raise Exception("Incorrect SNP frequencies {}, it should = {} !".format(snp_freq_sum, (row - n_c_row)))

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

ref_bases_sum = sum(ref_bases.values())
ref_bases_taken_sum = sum(ref_bases_taken.values())
if ref_bases_sum < ref_bases_taken_sum:
	raise Exception("Total nucleotides in VCF ref {} should >= totals in reference sequences {} !\n"
		"Please check ref_bases.json and ref_bases_taken.json".format(ref_bases_taken_sum, ref_bases_sum))

# count total constant site bases
constant_site_bases = ref_bases_sum - ref_bases_taken_sum
print("Total constant site bases =", constant_site_bases) 

# add constant site to patterns
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


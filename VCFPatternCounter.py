# Author: Walter Xie
# Last update: 04 Apr 2016

import vcf
import sys
import re
from Bio import SeqIO


# return a tuple of sample names on 2 copies 
def getSamples(vcf_reader):
	samples=() 
	#print (vcf_reader.samples)
	for sn in vcf_reader.samples:
		s1=sn + "_1"
		s2=sn + "_2"
		samples+=(s1,s2,)
	return samples

# return a string of one-site pattern in a site of all samples on 2 copies
def getOneSitePatttern(record, ref, alt): 
	print(record)
	site_patt=""
	for sample in record.samples:
		#sample = record.samples[1]
		gt=sample['GT']
		print(gt)
		if gt == ".":
			site_patt += ref * 2
		else:
			# ignore phase/unphase
			for gt2 in re.split('\||/', gt):
				gt2 = int(gt2)
				if gt2 > 0:
					site_patt += str(alt[gt2-1])
				elif gt2 == 0:
					site_patt += ref
				else:
					raise Exception("Invaild genotype, sample =", sample.sample, "GT =", gt)	
		
	if len(site_patt) != len(record.samples)*2:
		print(site_patt)
		raise Exception("One site pattern has incorrect number in", record.CHROM, "position", record.POS,
				"sample", sample.sample, len(site_patt), "!=", len(record.samples)*2)	
		
	return site_patt	


#vcf_f_n = sys.argv[1]
vcf_f_n="./adelie/allmodern.24.allancient.22.q20.vcf"
vcf_reader = vcf.Reader(open(vcf_f_n, 'r'))
print("Count patterns from", vcf_f_n)

rows_not_count_file = open("RowsNotCount.txt", "w")
patterns_file = open("Patterns.txt", "w")

samples=getSamples(vcf_reader)

patterns={}
row=0
for record in vcf_reader:
	#record = next(vcf_reader)
	#print (record)
	
	ref=record.REF
	# only pick up mutation, ignore insertion/deletion
	if len(ref)==1: 
		alt=record.ALT		
		site_patt=getOneSitePatttern(record, ref, alt)
		if site_patt in patterns:
			patterns[site_patt] += 1
		else:
			patterns[site_patt] = 1
				
		chrom=record.CHROM
		pos=record.POS
		
	else:
		print(record, file=rows_not_count_file)
	row += 1
	
print("Find", len(patterns), "patterns from total", row, "rows.")
print('\t'.join(samples), file=patterns_file)
for k, v in patterns.items():
	print('{}\t{}'.format(k, v), file=patterns_file)

patterns_file.close()
rows_not_count_file.close()


#ref_f_n = sys.argv[2]
ref_f_n="adelie.allscaffolds.fa"
reference_sequences = SeqIO.to_dict(SeqIO.parse(open(ref_f_n),'fasta'))
print("There are", len(reference_sequences), "sequences in the reference", ref_f_n)



def fib(n): 
	for ref_seq in reference_sequences:
		#ref_seq = next(reference_sequences)
		print(ref_seq.id)
		print(ref_seq.seq)
		print(len(ref_seq))






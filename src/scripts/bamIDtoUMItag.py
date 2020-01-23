#!/usr/bin/env python3

import pysam
import tempfile
import regex
from Bio import bgzf
import argparse

def barcode_to_tag(input_file, barcode, verbose):
	samfile = pysam.AlignmentFile(input_file, "rb")
	header = str(samfile.header)

	barcode_length = len(regex.findall('B', barcode))
	umi_length = len(regex.findall('U', barcode))
	barcode_pattern = '('
	if barcode_length > 0:
		barcode_pattern += '.BC:Z:[A-Z]{' + str(barcode_length) + '}'
	if umi_length > 0:
		barcode_pattern += '.RX:Z:[A-Z]{' + str(umi_length) + '}'
	barcode_pattern += ')'

	total = 0
	wrote = 0
	with tempfile.TemporaryFile() as tmp:
		with pysam.AlignmentFile(input_file, "rb") as infile:
			for read in infile:
				name_list = regex.split(barcode_pattern, read.query_name)
				read.query_name = ''.join([name_list[0], name_list[2]])
				tags = tuple(name_list[1].replace('.', ':Z:').split(':Z:'))[1:]
				read.tags = read.tags + [ tags[x:x + 2] for x in range(0, len(tags), 2) ]
				tmp.write((pysam.AlignedSegment.to_string(read)+'\n').encode('utf8'))
				total += 1
		tmp.seek(0)
		with bgzf.BgzfWriter(input_file, "wb") as outfile:
			outfile.write(header.encode('utf8'))
			for read in tmp:
				outfile.write(read)
				wrote += 1

	if verbose:
		print(total, "entries read from the input file.")
		print(wrote, "entries written to the output file.")


	return 0

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description="Generate the preseq input histogram", prefix_chars="-")
	parser.add_argument("-i", "--input", type=str, dest='infile', nargs=1, help="the input BAM file (not filtered by mapq)", required=True)
	parser.add_argument("-b", "--barcode", type=str, dest='barcode', nargs=1, help="the pattern of 'B' and 'U' characters identifying the barcode and/or UMI positions (e.g., BBBUUUBB)", required=True)
	parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true", default=False)
	args = parser.parse_args()

	input_file = args.infile[0]
	barcode = args.barcode[0]
	verbose = args.verbose

	barcode_to_tag(input_file, barcode, verbose)


#!/usr/bin/env python3

import sys
import gzip
import regex
import tempfile
from Bio.SeqIO.QualityIO import FastqGeneralIterator

import argparse

def move_barcode(input_file, output_file, skip5, barcode, verbose):
	total = 0
	wrote = 0

	barcode_length = len(barcode)
	barcode_list = [ i.start() for i in regex.finditer('B', barcode) ]
	umi_list = [ i.start() for i in regex.finditer('U', barcode) ]

	with tempfile.TemporaryFile() as tmp:
		with gzip.open(input_file[0], "rt") as in_handle:
			for title, seq, qual in FastqGeneralIterator(in_handle):
				seq = seq[skip5:]
				qual = qual[skip5:]
				tags = []
				if barcode_list:
					tags = tags + [":".join(["BC", "Z"] + [ "".join([ seq[i] for i in barcode_list ]) ])]
				if umi_list:
					tags = tags + [":".join(["RX", "Z"] + [ "".join([ seq[i] for i in umi_list ]) ])]
				title_list = title.strip().split(" ")[0].split(".")
				entry = "".join(['@' , (".".join(title_list[:-1] + tags + title_list[-1:])) , '\n' , seq[barcode_length:] , '\n+\n' + qual[barcode_length:] , '\n'])
				tmp.write(entry.encode('utf8'))
				total += 1
		
		tmp.seek(0)
		with gzip.open(output_file[0], "wt") as out_handle:
			for entry in tmp:
				out_handle.write(entry.decode('utf8'))
				wrote += 1

	if len(input_file) > 1:
		wrote_2 = 0
		with tempfile.TemporaryFile() as tmp:
			with gzip.open(output_file[0], "rt") as handle_1:
				with gzip.open(input_file[1], "rt") as handle_2:
					while True:
						line_1 = handle_1.readline()
						line_2 = handle_2.readline()
						if not (line_1 or line_2):
							break
						if line_1[0] == '@' or line_1[0] == '+':
							tmp.write(line_1.encode('utf8'))
						else:
							tmp.write(line_2.encode('utf8'))
			
			tmp.seek(0)
			with gzip.open(output_file[1], "wt") as out_handle:
				for entry in tmp:
					out_handle.write(entry.decode('utf8'))
					wrote_2 += 1

	if verbose:
		print(total, "entries read from the input file.")
		print(int(wrote/4), "entries written to the output file.")
		if len(input_file) > 1:
			print(int(wrote_2/4), "entries written to the output file 2.")

	return 0

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description="Pre-process the fastq file", prefix_chars="-")
	parser.add_argument("-i", "--input", type=str, dest='infile', nargs='+', help="the input fastq file (must be in .gz format)", required=True)
	parser.add_argument("-o", "--output", type=str, dest='outfile', nargs='+', help="the output fastq file (will be in .gz format)")
	parser.add_argument("-b", "--barcode", type=str, dest='barcode', nargs=1, help="the pattern of 'B' and 'U' characters identifying the barcode and/or UMI positions (e.g., BBBUUUBB)", required=True)
	parser.add_argument("-5", "--skip5", type=int, choices=list(range(0, 10)), default=[0], dest='skip5', nargs=1, help="the number of nucleotides to skip before the barcode/UMI (default: 0)")
	parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true", default=False)
	args = parser.parse_args()

	if len(args.infile) > 2:
		sys.exit('Wrong number of input file(s): it should be either 1 or 2!')
	if args.outfile:
		if len(args.outfile) != len(args.outfile):
			sys.exit('The number of input and output files must be the same!')
		if len(args.outfile) > 2:
			sys.exit('Wrong number of output file(s): it should be either 1 or 2!')

	input_file = args.infile
	if( args.outfile ):
		output_file = args.outfile
	else:
		output_file = input_file
	barcode = args.barcode[0]
	skip5 = args.skip5[0]
	verbose = args.verbose

	if verbose:
		print('Input file(s):', input_file)
		print('Output file(s):', output_file)
		print('Barcode pattern', str(barcode), "nt")
		print('The first', str(skip5), 'nt will be ignored')

	move_barcode(input_file, output_file, skip5, barcode, verbose)

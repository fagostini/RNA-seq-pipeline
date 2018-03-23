"""
An alternate method to pair fastq files
"""

import os
import sys

import argparse

def stream_fastq(fqfile):
	"""Read a fastq file and provide an iterable of the sequence ID, the
	full header, the sequence, and the quaity scores.

	Note that the sequence ID is the header up until the first space,
	while the header is the whole header.
	"""

	qin = open(fqfile, 'r')

	while True:
		header = qin.readline()
		if not header:
			break
		header = header.strip()
		seqidparts = header.split(' ')
		seqid = seqidparts[0]
		seq = qin.readline()
		seq = seq.strip()
		qualheader = qin.readline()
		qualscores = qin.readline()
		qualscores = qualscores.strip()
		header = header.replace('@', '', 1)
		yield seqid, header, seq, qualscores


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Pair fastq files, writing all the pairs to separate files and the unmapped reads to separate files")
	parser.add_argument('-l', help='Pair #1 reads file')
	parser.add_argument('-r', help='Pair #2 reads file')
	args = parser.parse_args()


	# read the first file into a data structure
	seqs = {}
	for (seqid, header, seq, qual) in stream_fastq(args.l):
		seqid = seqid.replace('.1', '').replace('/1', '')
		seqs[seqid] = [header, seq, qual]

	lp = open("{}.l.fq".format(args.l.replace('.1.fq', '')), 'w')
	rp = open("{}.r.fq".format(args.r.replace('.2.fq', '')), 'w')

	# read the first file into a data structure
	seen = set()
	for (seqid, header, seq, qual) in stream_fastq(args.r):
		seqid = seqid.replace('.2', '').replace('/2', '')
		seen.add(seqid)
		if seqid in seqs:
			lp.write("@" + seqs[seqid][0] + "\n" + seqs[seqid][1] + "\n+\n" + seqs[seqid][2] + "\n")
			rp.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")

	lp.close()
	rp.close()

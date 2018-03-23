"""
An alternate method to pair fastq files
"""

import os
import sys

import argparse

import gzip

def stream_fastq(fqfile):
	"""Read a fastq file and provide an iterable of the sequence ID, the
	full header, the sequence, and the quaity scores.

	Note that the sequence ID is the header up until the first space,
	while the header is the whole header.
	"""

	qin = gzip.open(fqfile, 'rb')

	while True:
		header = qin.readline()
		if not header:
			break
		header = header.strip()
		seqidparts = header.split(b' ')
		seqid = seqidparts[0]
		seq = qin.readline()
		seq = seq.strip()
		qualheader = qin.readline()
		qualscores = qin.readline()
		qualscores = qualscores.strip()
		header = header.replace(b'@', b'', 1)
		yield seqid, header, seq, qualscores


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Provided with two paired-end fastq files, the script matches the pairs overwriting the original files")
	parser.add_argument('-l', help='Pair #1 reads file')
	parser.add_argument('-r', help='Pair #2 reads file')
	args = parser.parse_args()


	# read the first file into a data structure
	seqs = {}
	for (seqid, header, seq, qual) in stream_fastq(args.l):
		seqid = seqid.replace(b'.1', b'').replace(b'/1', b'')
		seqs[seqid] = [header, seq, qual]

	lp = gzip.open("{}.l.fq.gz".format(args.l.replace('.fq.1.gz', '')), 'w')
	rp = gzip.open("{}.r.fq.gz".format(args.r.replace('.fq.2.gz', '')), 'w')

	# read the first file into a data structure
	seen = set()
	for (seqid, header, seq, qual) in stream_fastq(args.r):
		seqid = seqid.replace(b'.2', b'').replace(b'/2', b'')
		seen.add(seqid)
		if seqid in seqs:
			lp.write("@" + seqs[seqid][0] + "\n" + seqs[seqid][1] + "\n+\n" + seqs[seqid][2] + "\n")
			rp.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")

	lp.close()
	rp.close()

	os.rename("{}.l.fq.gz".format(args.l.replace('.fq.1.gz', '')), args.l)
	os.rename("{}.r.fq.gz".format(args.r.replace('.fq.2.gz', '')), args.r)

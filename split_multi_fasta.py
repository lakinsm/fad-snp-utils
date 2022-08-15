#!/usr/bin/env python3


import sys


def load_fasta(infile):
	with open(infile, 'r') as fasta_file:
		# Skip whitespace
		while True:
			line = fasta_file.readline()
			if line == "":
				return  # Empty file or premature end of file?
			if line[0] == ">":
				break
		while True:
			if line[0] != ">":
				raise ValueError("Records in FASTA should begin with '>'")
			header = line[1:].rstrip()
			all_lines = []
			line = fasta_file.readline()
			while True:
				if not line:
					break
				if line[0] == ">":
					break
				all_lines.append(line.rstrip())
				line = fasta_file.readline()
			yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
			if not line:
				return  # Stop Iteration


def output_fasta(F, outdir):
	for header, seq in F.items():
		with open(outdir + '/{}.fasta'.format(header.split('.')[0]), 'w') as out:
			out.write('>{}\n{}\n'.format(
				header.split('.')[0],
				seq
			))


if __name__ == '__main__':
	F = {k: v for k, v in load_fasta(sys.argv[1])}
	output_fasta(F, sys.argv[2])

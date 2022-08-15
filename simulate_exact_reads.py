#!/usr/bin/env python3

import sys
import os


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


def simulate_paired_reads(fastas, output_folder, samplename, insert_size=400, read_len=250, coverage=40):
	forward_out = open(output_folder + '/{}_S1_L001_R1_001.fastq'.format(samplename), 'w')
	reverse_out = open(output_folder + '/{}_S1_L001_R2_001.fastq'.format(samplename), 'w')
	transl_table = str.maketrans('ACGT', 'TGCA')
	counter = 1
	for header, seq in fastas.items():
		for _ in range(coverage):
			forward_out.write('@SimulatedExact:{} 1:N:0:2\n{}\n+\n{}\n'.format(
				counter,
				seq[:read_len],
				'F' * read_len
			))
			reverse_out.write('@SimulatedExact:{} 2:N:0:2\n{}\n+\n{}\n'.format(
				counter,
				seq[::-1][(insert_size - read_len):insert_size].translate(transl_table),
				'F' * read_len
			))
			counter += 1
	shift = insert_size // coverage

	for header, seq in fastas.items():
		for start in range(0, len(seq) - insert_size, shift):
			forward_out.write('@SimulatedExact:{} 1:N:0:2\n{}\n+\n{}\n'.format(
				counter,
				seq[start:(start + read_len)],
				'F' * read_len
			))
			reverse_out.write('@SimulatedExact:{} 2:N:0:2\n{}\n+\n{}\n'.format(
				counter,
				seq[::-1][(start + insert_size - read_len):(start + insert_size)].translate(transl_table),
				'F' * read_len
			))
			counter += 1
	for header, seq in fastas.items():
		for _ in range(coverage):
			forward_out.write('@SimulatedExact:{} 1:N:0:2\n{}\n+\n{}\n'.format(
				counter,
				seq[(len(seq) - insert_size):(len(seq) - insert_size + read_len)],
				'F' * read_len
			))
			reverse_out.write('@SimulatedExact:{} 2:N:0:2\n{}\n+\n{}\n'.format(
				counter,
				seq[::-1][(len(seq) - read_len):].translate(transl_table),
				'F' * read_len
			))
			counter += 1
	forward_out.close()
	reverse_out.close()


if __name__ == '__main__':
	F = {k: v for k, v in load_fasta(sys.argv[1])}
	output_dir = os.path.realpath(sys.argv[2])
	this_sample = sys.argv[3]
	simulate_paired_reads(F, output_dir, this_sample)

#!/usr/bin/env python

import argparse
import logging
import sys
from Bio.SeqUtils import nt_search
from Bio import SeqIO

choices = ['abi', 'ace', 'clustal', 'embl', 'fasta', 'fastq-sanger', 'fastq', 'fastq-solexa', 'fastq-illumina', 'genbank', 'gb', 'ig' \
			'imgt', 'nexus', 'phd', 'phylip', 'pir', 'seqxml', 'sff', 'stockholm', 'swiss', 'tab', 'qual', 'uniprot-xml']


def generate_fastas(fasta_sequences):
	"""
	Generator for extracting fasta sequences from a Bio.seqIO readable file.
	Returns the sequence and its ID.
	"""
	for f_seq in fasta_sequences:
		sequence, id = str(f_seq.seq), f_seq.id
		yield sequence, id


def find_pattern(query, pattern):
	"""
	Finds all occurrences of a pattern in the query sequence.
	Outputs sequence ID, start and end postion of the pattern.
	"""
	for (q,id) in query:
		find_match(q,pattern,id)
		find_match(q[::-1],pattern, id + 'R' )


def find_match(query,pattern,id):
	l = len(pattern)
	for match in nt_search(query, pattern):
		if type(match) == int:
			sys.stdout.write('{:30s} {:5d}  {:5d}'.format('>'+id,match,match+l)+'\n')
	return


if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-q', '--query' , required=True)
	parser.add_argument('-p', '--pattern' , required=True)
	parser.add_argument('-f', '--format', choices= choices)
	args = parser.parse_args()
	query = args.query
	pattern = args.pattern
	format = args.format
	
	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	logger.info('Searching  %s for %s. Forward and reverse strand search...'  %(query,pattern))
	
	find_pattern(generate_fastas(SeqIO.parse(open(query),format)),pattern)
	
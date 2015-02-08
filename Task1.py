
import argparse
import logging
import sys
from Bio.SeqUtils import nt_search
from Bio import SeqIO
import pandas as pd

"""------------------------------ Command-line Parser ----------------------------"""

parser = argparse.ArgumentParser()

parser.add_argument('-q', '--query' , required=True)
parser.add_argument('-p', '--pattern' , required=True)


args = parser.parse_args()
fq = args.query
p = args.pattern


"""------------------------------ Start the Logger -------------------------------"""


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info('Searching  %s for pattern:%s. Forward strand search...'  %(fq,p))


"""------------------------------ Functions---------------------------------------"""

def gf(fasta_sequences):
	for i in fasta_sequences:
		sequence, id = i.seq.tostring(), i.id
		yield sequence, id
		
def find_pattern(query, pattern):
	l = len(pattern)
	for (q,id) in query:
		for match in nt_search(q, pattern):
			if type(match) == int:
				sys.stdout.write('{:30s} {:5d}  {:5d}'.format('>'+id,match,match+l)+'\n')

find_pattern(gf(SeqIO.parse(open(fq),'fasta')),p)

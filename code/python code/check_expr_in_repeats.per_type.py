#!/usr/bin/env python3

"""
> check_expr_in_repeats.per_type.py <

Are dinoflagellate repeat elements differentially expressed with temperature?

(Code broadly follows check_meth_density_in_repeats.py.)
"""
import argparse
import collections
import csv
import gzip
import math
import pathlib
import re
import sys
import time

import numpy as np
import scipy.stats

import parse_fasta

parser = argparse.ArgumentParser(description="""
Are dinoflagellate repeat elements differentially expressed with temperature?""")

parser.add_argument('genome_fasta', metavar='fasta_file',
                    type=argparse.FileType('r'),
                    help='genome FASTA file.')
parser.add_argument('depth_tsv', metavar='tsv_filename',
                    type=argparse.FileType('r'), nargs='+',
                    help='Hisat2 depth .tsv filenames.')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

# hardcode repeats file location, as this is unlikely to change
REPEATS_FILE = '../2020-08-25-RepeatMasker/Mcavernosa_July2018.fasta.out.gz'

def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)

    seq = seq[::-1].translate(translation_table)

    return seq

# read sequences
genome_sequences = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')

if args.verbose:
    print (f'[{time.asctime()}] Finished reading {args.genome_fasta.name}...',
           file=sys.stderr)

# repeat context stores repeat types (LTR/LINE/SINE/...) from RepeatMasker.
#    0: not a repeat
#   1+: one of the repeat types
repeat_context = {}

# expr_context stores coverage values of underlying loci
#    0: not expressed
#   1+: expressed
expr_context = {}

# read in repeat regions from the repeats file. if repeats overlap, the
# subsequent ones in the file "win out"--for ease of programming, really
for s in genome_sequences:
    repeat_context[s] = np.zeros(len(genome_sequences[s]), np.int8)

repeat_type_to_int = {}
repeat_int = 0
with gzip.open(REPEATS_FILE, 'rt') as f:
    # skip first three lines
    f.readline()
    f.readline()
    f.readline()
    
    for line in f:
        line = re.split(r'\s+', line.strip())
        
        scaf = line[4]
        start_pos = int(line[5]) - 1
        end_pos = int(line[6])
        repeat_type = line[10].split('/')[0]  # i.e. "LTR/Gypsy" --> "LTR"
        if repeat_type not in repeat_type_to_int:
            repeat_int += 1
            repeat_type_to_int[repeat_type] = repeat_int
        
        repeat_context[scaf][start_pos:end_pos] = repeat_type_to_int[repeat_type]

if args.verbose:
    print (f'[{time.asctime()}] Finished parsing repeat context...',
           file=sys.stderr)

# read in meth/unmeth regions from Bismark cov(s), on a per-file basis
for dt in args.depth_tsv:
    # convert `dt` into a pathlib Path to make manipulation easier
    dt = pathlib.Path(dt.name)
    
    # reset the meth patterns
    for s in genome_sequences:
        expr_context[s] = np.zeros(len(genome_sequences[s]), np.int32)
    
    # deal with gzip-compressed tsvs
    if dt.name[-3:] == '.gz':
        tsv_reader = csv.reader(gzip.open(dt, 'rt'),
                                delimiter='\t')
    else:
        tsv_reader = csv.reader(open(dt), delimiter='\t')
    
    # define dict for coverages to be tallied
    all_combined = {x: 0 for x in range(len(repeat_type_to_int) + 1)}
    
    # read file, and set expr_context to the coverage value
    for row in tsv_reader:
        scaf = row[0]
        pos = int(row[1]) - 1   # convert to 0-based numbering
        coverage = int(row[2])
        
        # check which repeat context the position is from
        rc = repeat_context[scaf][pos]
        
        # increment coverage
        all_combined[rc] += coverage
    
    if args.verbose:
        print (f'[{time.asctime()}] {dt.name}: tallied everything up.',
               file=sys.stderr)
    
    # useful numbers
    sum_coverage = sum(all_combined.values())
    
    # print stuff out
    with open(dt.stem.split('.')[0] + '.expr.tsv', 'w') as f:
        print ('repeat_family', 'coverage', 'pct_coverage',
               sep='\t', file=f)
        
        # print a line for each repeat type
        for t in sorted(repeat_type_to_int):
            i = repeat_type_to_int[t]
            print (t, all_combined[i], all_combined[i] / sum_coverage * 100,
                   sep='\t', file=f)
        
        # print a line for not-repeats
        print ('not_repeats', all_combined[0], all_combined[0] / sum_coverage * 100,
               sep='\t', file=f)
    
    if args.verbose:
        print (f'[{time.asctime()}] {dt.name}: reporting done.',
               file=sys.stderr)

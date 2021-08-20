#!/usr/bin/env python3

"""
> check_meth_level_in_repeats.py <

If dinoflagellate CpGs/CHGs/CHHs are preferentially located within repeat
regions, are they also, on average, higher methylated?

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
If dinoflagellate CpGs/CHGs/CHHs are preferentially located within repeat
regions, are they also, on average, higher methylated?""")

parser.add_argument('genome_fasta', metavar='fasta_file',
                    type=argparse.FileType('r'),
                    help='genome FASTA file.')
parser.add_argument('bismark_cov', metavar='cov_filename',
                    type=argparse.FileType('r'), nargs='+',
                    help='Bismark .cov filename.')
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

def combine_triples(list_of_triples):
    """
    When given a bunch of triples, merge them together.
    """
    # in case the list is a generator, convert it to a list
    list_of_triples = [x for x in list_of_triples]
    
    combined_n = sum(x[0] for x in list_of_triples)
    combined_x = sum(x[1] for x in list_of_triples)
    combined_x2 = sum(x[2] for x in list_of_triples)
    
    return [combined_n, combined_x, combined_x2]

def compute_mean(n_sumx_sumx2):
    """
    Function takes in a single list with three numbers, and computes mean.
    """
    n = n_sumx_sumx2[0]
    sumx = n_sumx_sumx2[1]
    
    if n < 1: return 'NA'
    
    return sumx / n

def compute_sd(n_sumx_sumx2):
    """
    Function takes in a single list with three numbers, and computes 
    standard deviation (SD).
    """
    n = n_sumx_sumx2[0]
    sumx = n_sumx_sumx2[1]
    sumx2 = n_sumx_sumx2[2]
    
    if n < 2: return 'NA'
    
    samp_var = (sumx2 - (sumx ** 2) / n) / (n - 1)
    samp_sd = math.sqrt(samp_var)
    return samp_sd

def compute_sem(n_sumx_sumx2):
    """
    Function takes in a single list with three numbers, and computes 
    standard error of the mean (SEM).
    """
    n = n_sumx_sumx2[0]
    sumx = n_sumx_sumx2[1]
    sumx2 = n_sumx_sumx2[2]
    
    if n < 2: return 'NA'
    
    samp_var = (sumx2 - (sumx ** 2) / n) / (n - 1)
    samp_sd = math.sqrt(samp_var)
    return samp_sd / math.sqrt(n)

def perform_t_test(first_triple, second_triple):
    """
    Perform a t-test on `first_triple` vs. `second_triple`.
    """
    n1 = first_triple[0]
    n2 = second_triple[0]
    if n1 < 2 or n2 < 2: return 'NA'
    
    mean1 = first_triple[1] / n1
    std1 = compute_sd(first_triple)
    mean2 = second_triple[1] / n2
    std2 = compute_sd(second_triple)
    
    return scipy.stats.ttest_ind_from_stats(mean1, std1, n1,
                                            mean2, std2, n2, equal_var=False)[1]

def subtract_triples(first_triple, second_triple):
    """
    Carry out `first_triple` - `second_triple`.
    """
    return [first_triple[0] - second_triple[0],
            first_triple[1] - second_triple[1],
            first_triple[2] - second_triple[2]]

# read sequences
genome_sequences = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')

if args.verbose:
    print (f'[{time.asctime()}] Finished reading {args.genome_fasta.name}...',
           file=sys.stderr)

# genomic context is denoted in a NumPy array of int8s as follows:
#   0: not a cytosine on Watson/Crick strand
#   1: cytosine in CpG context on Watson/Crick strand
#   2: cytosine in CHG context on Watson/Crick strand
#   3: cytosine in CHH context on Watson/Crick strand
genomic_context = {}

# repeat context stores repeat types (LTR/LINE/SINE/...) from RepeatMasker.
#    0: not a repeat
#   1+: one of the repeat types
repeat_context = {}

# meth_context stores beta value of underlying loci
#          0: not methylated
#   (0, 100]: methylated
meth_context = {}

# read genomic context
for s in genome_sequences:
    genomic_context[s] = np.zeros(len(genome_sequences[s]), np.int8)
    repeat_context[s] = np.zeros(len(genome_sequences[s]), np.int8)
    
    scaf_watson = genome_sequences[s]
    scaf_crick = reverse_complement(genome_sequences[s])
    
    # get CpG
    for x in re.finditer('CG', scaf_watson):
        genomic_context[s][x.span()[0]] = 1
    
    for x in re.finditer('CG', scaf_crick):
        genomic_context[s][-x.span()[0] - 1] = 1
    
    # get CHG
    for x in re.finditer('C[A|C|T]G', scaf_watson):
        genomic_context[s][x.span()[0]] = 2
    
    for x in re.finditer('C[A|C|T]G', scaf_crick):
        genomic_context[s][-x.span()[0] - 1] = 2
    
    # get CHH--note the lookahead assertion: by default, Python only does non-
    # overlapping matches, which causes stuff like CCCC to report only 1 CHH.
    # lookahead assertion would correctly report 2 CHHs in CCCC.
    for x in re.finditer('(?=C[A|C|T][A|C|T])', scaf_watson):
        genomic_context[s][x.span()[0]] = 3
    
    for x in re.finditer('(?=C[A|C|T][A|C|T])', scaf_crick):
        genomic_context[s][-x.span()[0] - 1] = 3

if args.verbose:
    print (f'[{time.asctime()}] Finished parsing genomic context...',
           file=sys.stderr)

# read in repeat regions from the repeats file. if repeats overlap, the
# subsequent ones in the file "win out"--for ease of programming, really
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
for bc in args.bismark_cov:
    # convert bc into a pathlib Path to make manipulation easier
    bc = pathlib.Path(bc.name)
    
    # reset the meth patterns
    for s in genome_sequences:
        meth_context[s] = np.zeros(len(genome_sequences[s]), np.float32)
    
    # deal with gzip-compressed covs
    if bc.name[-3:] == '.gz':
        tsv_reader = csv.reader(gzip.open(bc, 'rt'),
                                delimiter='\t')
    else:
        tsv_reader = csv.reader(open(bc), delimiter='\t')
    
    # read file, and set meth_context of pos to 1 if it's methylated
    for row in tsv_reader:
        scaf = row[0]
        pos = int(row[1]) - 1   # convert to 0-based numbering
        meth_context[scaf][pos] = float(row[3])
    
    # variance/sd can be calculated without storing the entire list of
    # numbers, but three pieces of information required: n, sum(x), sum(x ** 2)
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    all_combined = {(x, y): [0, 0, 0] for x in range(4) 
                                      for y in range(len(repeat_type_to_int) + 1)}
    
    # walk through genome and increment the three numbers respectively
    for scaf in genomic_context:
        for pos in range(len(genomic_context[scaf])):
            # skip non-methylated positions
            if not meth_context[scaf][pos]: continue
            
            # define tuple containing all three contexts
            tupl = (genomic_context[scaf][pos], repeat_context[scaf][pos])
            
            # increment n
            all_combined[tupl][0] += 1
            
            # increment sum(x)
            all_combined[tupl][1] += meth_context[scaf][pos]
            
            # increment sum(x ** 2)
            all_combined[tupl][2] += meth_context[scaf][pos] ** 2
    
    if args.verbose:
        print (f'[{time.asctime()}] {bc.name}: tallied everything up.',
               file=sys.stderr)
    
    # useful numbers
    mcpg_all = combine_triples(all_combined[(x, y)] for x, y in all_combined if x == 1)
    mchg_all = combine_triples(all_combined[(x, y)] for x, y in all_combined if x == 2)
    mchh_all = combine_triples(all_combined[(x, y)] for x, y in all_combined if x == 3)
    mc_all = combine_triples([mcpg_all, mchg_all, mchh_all])
    
    # print stuff out
    with open(bc.stem.replace('.cov', '') + '.meth_level.tsv', 'w') as f:
        print ('repeat_family',
               'C_meth_in_repeat', 'n', 'sem',
               'C_meth_not_in_repeat', 'n', 'sem',
               'ttest_C',
               'CpG_meth_in_repeat', 'n', 'sem',
               'CpG_meth_not_in_repeat', 'n', 'sem',
               'ttest_CpG',
               'CHG_meth_in_repeat', 'n', 'sem',
               'CHG_meth_not_in_repeat', 'n', 'sem',
               'ttest_CHG',
               'CHH_meth_in_repeat', 'n', 'sem',
               'CHH_meth_not_in_repeat', 'n', 'sem',
               'ttest_CHH',
               sep='\t', file=f)
        
        # print a line for each repeat type
        for t in sorted(repeat_type_to_int):
            i = repeat_type_to_int[t]
            
            mcpg_in_rep = all_combined[(1, i)]
            mchg_in_rep = all_combined[(2, i)]
            mchh_in_rep = all_combined[(3, i)]
            mc_in_rep = combine_triples([mcpg_in_rep, mchg_in_rep, mchh_in_rep])
            
            mcpg_not_in_rep = subtract_triples(mcpg_all, mcpg_in_rep)
            mchg_not_in_rep = subtract_triples(mchg_all, mchg_in_rep)
            mchh_not_in_rep = subtract_triples(mchh_all, mchh_in_rep)
            mc_not_in_rep = combine_triples([mcpg_not_in_rep, mchg_not_in_rep, mchh_not_in_rep])
            
            print (t,
                   compute_mean(mc_in_rep), mc_in_rep[0], compute_sem(mc_in_rep),
                   compute_mean(mc_not_in_rep), mc_not_in_rep[0], compute_sem(mc_not_in_rep),
                   perform_t_test(mc_in_rep, mc_not_in_rep),
                   compute_mean(mcpg_in_rep), mcpg_in_rep[0], compute_sem(mcpg_in_rep),
                   compute_mean(mcpg_not_in_rep), mcpg_not_in_rep[0], compute_sem(mcpg_not_in_rep),
                   perform_t_test(mcpg_in_rep, mcpg_not_in_rep),
                   compute_mean(mchg_in_rep), mchg_in_rep[0], compute_sem(mchg_in_rep),
                   compute_mean(mchg_not_in_rep), mchg_not_in_rep[0], compute_sem(mchg_not_in_rep),
                   perform_t_test(mchg_in_rep, mchg_not_in_rep),
                   compute_mean(mchh_in_rep), mchh_in_rep[0], compute_sem(mchh_in_rep),
                   compute_mean(mchh_not_in_rep), mchh_not_in_rep[0], compute_sem(mchh_not_in_rep),
                   perform_t_test(mchh_in_rep, mchh_not_in_rep),
                   sep='\t', file=f)
        
        # and print a final line that sums everything up
        mcpg_in_rep = combine_triples(all_combined[(x, y)] for x, y in all_combined 
                                      if x == 1 and y > 0)
        mchg_in_rep = combine_triples(all_combined[(x, y)] for x, y in all_combined 
                                      if x == 2 and y > 0)
        mchh_in_rep = combine_triples(all_combined[(x, y)] for x, y in all_combined 
                                      if x == 3 and y > 0)
        mc_in_rep = combine_triples([mcpg_in_rep, mchg_in_rep, mchh_in_rep])
        
        mcpg_not_in_rep = subtract_triples(mcpg_all, mcpg_in_rep)
        mchg_not_in_rep = subtract_triples(mchg_all, mchg_in_rep)
        mchh_not_in_rep = subtract_triples(mchh_all, mchh_in_rep)
        mc_not_in_rep = combine_triples([mcpg_not_in_rep, mchg_not_in_rep, mchh_not_in_rep])
        
        print ('all_repeats',
               compute_mean(mc_in_rep), mc_in_rep[0], compute_sem(mc_in_rep),
               compute_mean(mc_not_in_rep), mc_not_in_rep[0], compute_sem(mc_not_in_rep),
               perform_t_test(mc_in_rep, mc_not_in_rep),
               compute_mean(mcpg_in_rep), mcpg_in_rep[0], compute_sem(mcpg_in_rep),
               compute_mean(mcpg_not_in_rep), mcpg_not_in_rep[0], compute_sem(mcpg_not_in_rep),
               perform_t_test(mcpg_in_rep, mcpg_not_in_rep),
               compute_mean(mchg_in_rep), mchg_in_rep[0], compute_sem(mchg_in_rep),
               compute_mean(mchg_not_in_rep), mchg_not_in_rep[0], compute_sem(mchg_not_in_rep),
               perform_t_test(mchg_in_rep, mchg_not_in_rep),
               compute_mean(mchh_in_rep), mchh_in_rep[0], compute_sem(mchh_in_rep),
               compute_mean(mchh_not_in_rep), mchh_not_in_rep[0], compute_sem(mchh_not_in_rep),
               perform_t_test(mchh_in_rep, mchh_not_in_rep),
               sep='\t', file=f)
    
    if args.verbose:
        print (f'[{time.asctime()}] {bc.name}: reporting done.',
               file=sys.stderr)

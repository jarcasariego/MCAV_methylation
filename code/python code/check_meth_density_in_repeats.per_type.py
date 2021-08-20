#!/usr/bin/env python3

"""
> check_meth_density_in_repeats.py <

Are dinoflagellate CpGs/CHGs/CHHs preferentially located within repeat regions,
a la plants?
"""
import argparse
import collections
import csv
import gzip
import pathlib
import re
import sys
import time

import numpy as np
import scipy.stats

import parse_fasta

parser = argparse.ArgumentParser(description="""
Are dinoflagellate CpGs/CHGs/CHHs preferentially located within repeat regions,
a la plants?""")

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
# it's another NumPy array of int8s, but starting from 1
repeat_context = {}

# meth_context stores whether the underlying loci is methylated
#   0: not methylated
#   1: methylated
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
# exclude ARTEFACT, Retroposon and scRNA repeat regions 
# individual repeat element types
excluded_repeats = ['ARTEFACT', 'Retroposon', 'scRNA']
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
        if repeat_type not in excluded_repeats:
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
        meth_context[s] = np.zeros(len(genome_sequences[s]), np.int8)
    
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
        meth_context[scaf][pos] = 1
    
    # tally and print stuff out, on a per-file basis
    all_combined = {(x, y, z): 0 for x in range(4) 
                                 for y in range(len(repeat_type_to_int) + 1) 
                                 for z in range(2)}
    for scaf in genomic_context:
        for pos in range(len(genomic_context[scaf])):
            # append tuples of all three contexts
            tupl = (genomic_context[scaf][pos], repeat_context[scaf][pos],
                    meth_context[scaf][pos])
            all_combined[tupl] += 1
    
    if args.verbose:
        print (f'[{time.asctime()}] {bc.name}: tallied everything up.',
               file=sys.stderr)
    
    # useful numbers
    n_cpg = sum(all_combined[(x, y, z)] for x, y, z in all_combined if x == 1)
    n_chg = sum(all_combined[(x, y, z)] for x, y, z in all_combined if x == 2)
    n_chh = sum(all_combined[(x, y, z)] for x, y, z in all_combined if x == 3)
    n_c = n_cpg + n_chg + n_chh
    
    n_meth_cpg = sum(all_combined[(x, y, z)] for x, y, z in all_combined if x == 1 and z == 1)
    n_meth_chg = sum(all_combined[(x, y, z)] for x, y, z in all_combined if x == 2 and z == 1)
    n_meth_chh = sum(all_combined[(x, y, z)] for x, y, z in all_combined if x == 3 and z == 1)
    n_meth_c = n_meth_cpg + n_meth_chg + n_meth_chh
    
    n_unmeth_cpg = sum(all_combined[(x, y, z)] for x, y, z in all_combined if x == 1 and z == 0)
    n_unmeth_chg = sum(all_combined[(x, y, z)] for x, y, z in all_combined if x == 2 and z == 0)
    n_unmeth_chh = sum(all_combined[(x, y, z)] for x, y, z in all_combined if x == 3 and z == 0)
    n_unmeth_c = n_unmeth_cpg + n_unmeth_chg + n_unmeth_chh
    
    # print stuff out
    with open(bc.stem.replace('.cov', '') + '.meth_density.tsv', 'w') as f:
        print ('repeat_family', 'meth_CpG', 'meth_CHG', 'meth_CHH',
               'unmeth_CpG', 'unmeth_CHG', 'unmeth_CHH', 
               'within_meth_C', 'outside_meth_C',               # }
               'odds_ratio_meth_C', 'fisher_p_meth_C',          # }
               'within_meth_CpG', 'outside_meth_CpG',           # }
               'odds_ratio_meth_CpG', 'fisher_p_meth_CpG',      # } stats stuff
               'within_meth_CHG', 'outside_meth_CHG',           # }
               'odds_ratio_meth_CHG', 'fisher_p_meth_CHG',      # }
               'within_meth_CHH', 'outside_meth_CHH',           # }
               'odds_ratio_meth_CHH', 'fisher_p_meth_CHH',      # }
               sep='\t', file=f)
        
        # print a line for each repeat type
        for t in sorted(repeat_type_to_int):
            i = repeat_type_to_int[t]
            
            mcpg_in_rep = all_combined[(1, i, 1)]
            mchg_in_rep = all_combined[(2, i, 1)]
            mchh_in_rep = all_combined[(3, i, 1)]
            mc_in_rep = mcpg_in_rep + mchg_in_rep + mchh_in_rep
            
            mcpg_not_in_rep = n_meth_cpg - mcpg_in_rep
            mchg_not_in_rep = n_meth_chg - mchg_in_rep
            mchh_not_in_rep = n_meth_chh - mchh_in_rep
            mc_not_in_rep = n_meth_c - mc_in_rep
            
            umcpg_in_rep = all_combined[(1, i, 0)]
            umchg_in_rep = all_combined[(2, i, 0)]
            umchh_in_rep = all_combined[(3, i, 0)]
            umc_in_rep = umcpg_in_rep + umchg_in_rep + umchh_in_rep
            
            umcpg_not_in_rep = n_unmeth_cpg - umcpg_in_rep
            umchg_not_in_rep = n_unmeth_chg - umchg_in_rep
            umchh_not_in_rep = n_unmeth_chh - umchh_in_rep
            umc_not_in_rep = n_unmeth_c - umc_in_rep
            
            print (t, mcpg_in_rep, mchg_in_rep, mchh_in_rep,
                   umcpg_in_rep, umchg_in_rep, umchh_in_rep,
                   mc_in_rep / n_meth_c, mc_not_in_rep / n_meth_c,
                   (mc_in_rep / umc_in_rep) / (mc_not_in_rep / umc_not_in_rep),
                   scipy.stats.fisher_exact([
                       [mc_in_rep, umc_in_rep], [mc_not_in_rep, umc_not_in_rep]])[1],
                   mcpg_in_rep / n_meth_cpg, mcpg_not_in_rep / n_meth_cpg,
                   (mcpg_in_rep / umcpg_in_rep) / (mcpg_not_in_rep / umcpg_not_in_rep),
                   scipy.stats.fisher_exact([
                       [mcpg_in_rep, umcpg_in_rep], [mcpg_not_in_rep, umcpg_not_in_rep]])[1],
                   mchg_in_rep / n_meth_chg, mchg_not_in_rep / n_meth_chg,
                   (mchg_in_rep / umchg_in_rep) / (mchg_not_in_rep / umchg_not_in_rep),
                   scipy.stats.fisher_exact([
                       [mchg_in_rep, umchg_in_rep], [mchg_not_in_rep, umchg_not_in_rep]])[1],
                   mchh_in_rep / n_meth_chh, mchh_not_in_rep / n_meth_chh,
                   (mchh_in_rep / umchh_in_rep) / (mchh_not_in_rep / umchh_not_in_rep),
                   scipy.stats.fisher_exact([
                       [mchh_in_rep, umchh_in_rep], [mchh_not_in_rep, umchh_not_in_rep]])[1],
                   sep='\t', file=f)
        
        # and print a final line that sums everything up
        mcpg_in_rep = sum(all_combined[(x, y, z)] for x, y, z in all_combined 
                          if x == 1 and y > 0 and z == 1)
        mchg_in_rep = sum(all_combined[(x, y, z)] for x, y, z in all_combined
                          if x == 2 and y > 0 and z == 1)
        mchh_in_rep = sum(all_combined[(x, y, z)] for x, y, z in all_combined
                          if x == 3 and y > 0 and z == 1)
        mc_in_rep = mcpg_in_rep + mchg_in_rep + mchh_in_rep
        
        mcpg_not_in_rep = n_meth_cpg - mcpg_in_rep
        mchg_not_in_rep = n_meth_chg - mchg_in_rep
        mchh_not_in_rep = n_meth_chh - mchh_in_rep
        mc_not_in_rep = n_meth_c - mc_in_rep
        
        umcpg_in_rep = sum(all_combined[(x, y, z)] for x, y, z in all_combined
                           if x == 1 and y > 0 and z == 0)
        umchg_in_rep = sum(all_combined[(x, y, z)] for x, y, z in all_combined
                           if x == 2 and y > 0 and z == 0)
        umchh_in_rep = sum(all_combined[(x, y, z)] for x, y, z in all_combined
                           if x == 3 and y > 0 and z == 0)
        umc_in_rep = umcpg_in_rep + umchg_in_rep + umchh_in_rep
        
        umcpg_not_in_rep = n_unmeth_cpg - umcpg_in_rep
        umchg_not_in_rep = n_unmeth_chg - umchg_in_rep
        umchh_not_in_rep = n_unmeth_chh - umchh_in_rep
        umc_not_in_rep = n_unmeth_c - umc_in_rep
        
        print ('sum_repeats', mcpg_in_rep, mchg_in_rep, mchh_in_rep,
               umcpg_in_rep, umchg_in_rep, umchh_in_rep,
               mc_in_rep / n_meth_c, mc_not_in_rep / n_meth_c,
               (mc_in_rep / umc_in_rep) / (mc_not_in_rep / umc_not_in_rep),
               scipy.stats.fisher_exact([
                   [mc_in_rep, umc_in_rep], [mc_not_in_rep, umc_not_in_rep]])[1],
               mcpg_in_rep / n_meth_cpg, mcpg_not_in_rep / n_meth_cpg,
               (mcpg_in_rep / umcpg_in_rep) / (mcpg_not_in_rep / umcpg_not_in_rep),
               scipy.stats.fisher_exact([
                   [mcpg_in_rep, umcpg_in_rep], [mcpg_not_in_rep, umcpg_not_in_rep]])[1],
               mchg_in_rep / n_meth_chg, mchg_not_in_rep / n_meth_chg,
               (mchg_in_rep / umchg_in_rep) / (mchg_not_in_rep / umchg_not_in_rep),
               scipy.stats.fisher_exact([
                   [mchg_in_rep, umchg_in_rep], [mchg_not_in_rep, umchg_not_in_rep]])[1],
               mchh_in_rep / n_meth_chh, mchh_not_in_rep / n_meth_chh,
               (mchh_in_rep / umchh_in_rep) / (mchh_not_in_rep / umchh_not_in_rep),
               scipy.stats.fisher_exact([
                   [mchh_in_rep, umchh_in_rep], [mchh_not_in_rep, umchh_not_in_rep]])[1],
               sep='\t', file=f)
        
        # and print a final line for not-repeats
        print ('not_repeats', 
               all_combined[(1, 0, 1)],
               all_combined[(2, 0, 1)],
               all_combined[(3, 0, 1)],
               all_combined[(1, 0, 0)],
               all_combined[(2, 0, 0)],
               all_combined[(3, 0, 0)], '\t' * 16, sep='\t', file=f)
    
    if args.verbose:
        print (f'[{time.asctime()}] {bc.name}: reporting done.',
               file=sys.stderr)

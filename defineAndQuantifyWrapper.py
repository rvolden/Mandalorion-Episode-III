#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', type=str,
                    help='Tab delimited file that specifies where minimap, blat, \
                          emtrey, and racon executables are')
parser.add_argument('-p', '--path', type=str, help='Directory to put output files into')
parser.add_argument('-m', '--score_matrix', type=str, default='./NUC.4.4.mat',
                    help='Substitution matrix (default ./NUC.4.4.mat)')
parser.add_argument('-u', '--upstream_buffer', type=str, default='10',
                    help='Defines leniency window for TSS determination (default 10)')
parser.add_argument('-d', '--downstream_buffer', type=str, default='50',
                    help='Defines leniency window for polyA determination (default 50)')
parser.add_argument('-s', '--subsample_consensus', type=str, default='500',
                    help='Defines how many random subreads are used to make isoforms (default 500)')
parser.add_argument('-g', '--genome_annotation', type=str,
                    help='Genome annotation file (gtf)')
parser.add_argument('-G', '--genome_sequence', type=str,
                    help='Genome file (fasta)')
parser.add_argument('-r', '--minimum_ratio', type=str, default='0.05',
                    help='Proportion of reads that align to a locus required for an \
                          isoform (default 0.05)')
parser.add_argument('-i', '--minimum_internal_ratio', type=str, default='0.125')
parser.add_argument('-R', '--minimum_reads', type=str, default='5',
                    help='Minimum number of reads to make an isoform (default 5)')
parser.add_argument('-a', '--adapter_file', type=str,
                    help='Fasta file with 5prime and 3prime adapters')
parser.add_argument('-f', '--R2C2_Consensus_reads', type=str,
                    help='Fasta file with R2C2 consensus reads')
parser.add_argument('-b', '--R2C2_subreads', type=str,
                    help='Fastq file with R2C2 subreads')
parser.add_argument('-O', '--overhangs', type=str, default='0,40,0,40',
                    help='Defines bounds for unaligned bases on ends. \
                          Format: min5prime,max5prime,min3prime,max3prime (default 0,40,0,40))')
parser.add_argument('-t', '--minimap2_threads', type=str, default='4',
                    help='Number of threads to use when running minimap (default 4)')
parser.add_argument('-e', '--ends', type=str, default='ATGGG,AAAAA',
                    help='Ends of your sequences. Defaults to Smartseq ends.\
                          Format: 5prime,3prime')

args = parser.parse_args()

config_file = args.config_file
path = args.path + '/'	         #path where you want your output files to go
matrix = args.score_matrix
upstream_buffer = args.upstream_buffer
downstream_buffer = args.downstream_buffer
subsample_consensus = args.subsample_consensus
genome_annotation = args.genome_annotation
genome_sequence = args.genome_sequence
adapter = args.adapter_file
minimum_ratio = args.minimum_ratio
minimum_internal_ratio = args.minimum_internal_ratio
minimum_reads = args.minimum_reads
fasta_file = args.R2C2_Consensus_reads
subreads = args.R2C2_subreads
overhangs = args.overhangs
minimap2_threads = args.minimap2_threads
ends = args.ends

def configReader(configIn):
    '''Parses the config file.'''
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]
    # should have minimap, racon, consensus, blat, and emtrey
    possible = set(['minimap2', 'consensus', 'racon', 'blat', 'emtrey'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
    # check for missing programs
    # if missing, default to path
    for missing in possible-inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs

progs = configReader(config_file)
minimap2 = progs['minimap2']
racon = progs['racon']
consensus = progs['consensus']
emtrey = progs['emtrey']
consensus = 'python3 ' + consensus

print('Aligning reads')
sam_file = fasta_file + '.sam'
psl_file = fasta_file + '.psl'
clean_psl_file = fasta_file + '.clean.psl'
os.system('%s --secondary=no -ax splice -t %s %s %s > %s ' %(minimap2, minimap2_threads, genome_sequence, fasta_file, sam_file))
os.system('%s -i %s > %s ' %(emtrey, sam_file, psl_file))
os.system('%s %s %s ' %('python3 clean_psl.py', psl_file, clean_psl_file))
print('Finding Splice sites')
os.system('python3 spliceSites.py %s %s %s %s %s %s' %(clean_psl_file, path, '0.05', genome_annotation, 'g', sam_file))
print('Identifying Isoforms')
os.system('python3 defineAndQuantifyIsoforms.py %s %s %s %s %s %s' %(clean_psl_file, path, downstream_buffer, upstream_buffer, subreads, fasta_file)) # This script sort raw reads into isoform bins. The two number variables determine the window around TSS and TES in which read ends can fall and still be matched to the site.
os.system('python3 createConsensi.py -p %s -s %s -c %s -m %s' %(path, subsample_consensus, config_file, matrix))
os.system('python3 filterIsoforms.py -p %s -i %s -r %s -R %s -n %s -a %s -G %s -c %s -O %s -t %s -e %s' %(path, path+'/Isoform_Consensi.fasta', minimum_ratio, minimum_reads, minimum_internal_ratio, adapter, genome_sequence, config_file, overhangs, minimap2_threads, ends))

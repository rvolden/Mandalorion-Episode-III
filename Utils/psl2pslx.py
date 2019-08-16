import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--psl_file', type=str)
parser.add_argument('-g', '--genome_file', type=str)
parser.add_argument('-r', '--read_file', type=str)
parser.add_argument('-x', '--pslx_file', type=str)

args = parser.parse_args()
psl_file = args.psl_file
genome_file = args.genome_file
read_file = args.read_file
pslx_file= args.pslx_file

def reverse_complement(sequence):
    '''Returns the reverse complement of a sequence'''
    bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
    return ''.join([bases[x] for x in list(sequence)])[::-1]

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            readDict[line[1:]] = ''
            lastHead = line[1:]
        else:
            readDict[lastHead] += line.upper()
    return readDict

def parse_contigs(psl_file, pslx_file, genome, reads):
    out = open(pslx_file, 'w')
    for line in open(psl_file):
        a = line.strip().split('\t')
        chromosome = a[13]
        name = a[9]
        direction = a[8]
        mismatches = 0
        genomic = genome[chromosome]
        read = reads[name]
        if direction == '-':
            read = reverse_complement(read)
        start = int(a[15])
        blocksizes = np.array(a[18].split(',')[:-1], dtype=int)
        blockstarts = np.array(a[20].split(',')[:-1], dtype=int)
        readstarts = np.array(a[19].split(',')[:-1], dtype=int)
        qseqs, tseqs = '', ''
        for index in range(0,len(blocksizes),1):
            blocksize = blocksizes[index]
            blockstart = blockstarts[index]
            readstart = readstarts[index]

            read_block_sequence = read[readstart : readstart+blocksize]
            genome_block_sequence = genomic[blockstart : blockstart+blocksize]
            for i in range(blocksize):
                readBase = read_block_sequence[i]
                genomeBase = genome_block_sequence[i]
                if readBase != genomeBase:
                    mismatches += 1
            qseqs += read_block_sequence + ','
            tseqs += genome_block_sequence + ','

        a[1] = str(mismatches)
        new_line = ('\t').join(a)
        new_line += '\t' + qseqs + '\t' + tseqs
        out.write(new_line + '\n')

def main():
    print('reading genome')
    genome = read_fasta(genome_file)
    print('reading reads')
    reads = read_fasta(read_file)
    print('adding sequence info')
    parse_contigs(psl_file, pslx_file, genome, reads)

main()

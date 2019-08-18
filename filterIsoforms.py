#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import numpy as np
import argparse

def argParser():
    parser = argparse.ArgumentParser(description = 'Filters Isoforms \
                                                    generated by Mandalorion',
                                     add_help = True,
                                     prefix_chars = '-')

    parser.add_argument('--path', '-p', type=str, action='store', default=os.getcwd(),
                        help='Directory where all the files are/where they will end up.\
                              Defaults to your current directory.')
    parser.add_argument('--infile', '-i', type=str, action='store')
    parser.add_argument('--config', '-c', type=str, action='store', default='',
                        help='If you want to use a config file to specify paths to\
                              programs, specify them here. Use for poa, racon, gonk,\
                              blat, and minimap2 if they are not in your path.')
    parser.add_argument('--internal_ratio', '-n', type=float, action='store')
    parser.add_argument('-r', '--minimum_ratio', type=float)
    parser.add_argument('-R', '--minimum_reads', type=float)
    parser.add_argument('-G', '--genome_sequence', type=str)
    parser.add_argument('-a', '--adapter', type=str)
    parser.add_argument('-O', '--overhangs', type=str)
    parser.add_argument('-t', '--minimap2_threads', type=str)
    parser.add_argument('-e', '--ends', type=str, default='ATGGG,AAAAA',
                        help='Ends of your sequences. Defaults to Smartseq ends.\
                              Format: 5prime,3prime')

    return vars(parser.parse_args())

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

args = argParser()
path = args['path']
infile = args['infile']
minimum_ratio = args['minimum_ratio']
minimum_reads = args['minimum_reads']

internal_ratio = args['internal_ratio']

genome = args['genome_sequence']
configIn = args['config']
adapter = args['adapter']
overhangs = np.array(args['overhangs'].split(','), dtype=int)
minimap2_threads = args['minimap2_threads']
ends = args['ends']

if args['config']:
    progs = configReader(args['config'])
    minimap2 = progs['minimap2']
    racon = progs['racon']
    consensus = progs['consensus']
    emtrey = progs['emtrey']
    blat = progs['blat']
else:
    minimap2, racon, emtrey, blat = 'minimap2', 'racon', 'emtrey', 'blat'
    consensus = 'consensus.py'

consensus = 'python3 ' + consensus

out2 = open(path + '/Isoform_Consensi_filtered.fasta', 'w')
out3 = open(path + '/Isoform_Consensi_filtered.aligned.out.clean.psl', 'w')

def reverse_complement(sequence):
    '''Returns the reverse complement of a sequence'''
    bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
    return ''.join([bases[x] for x in list(sequence)])[::-1]

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    tempSeqs, headers, sequences = [], [], []
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            headers.append(line.split()[0][1:])
        # covers the case where the file ends while reading sequences
        if line.startswith('>'):
            sequences.append(''.join(tempSeqs).upper())
            tempSeqs = []
        else:
            tempSeqs.append(line)
    sequences.append(''.join(tempSeqs).upper())
    sequences = sequences[1:]
    for i in range(len(headers)):
        if headers[i].split('_')[-1] != '0' and sequences[i]:
            readDict[headers[i]] = sequences[i]
    return readDict

def get_count(isoform_list, chromosome, psl_dict):
    count = {}
    for isoform, info in psl_dict.items():
        coordinates = info[0]
        direction = info[1]
        number = int(isoform.split('_')[-1])
        counter = 0
        start, end = coordinates[0], coordinates[-1]
        for base in np.arange(round(start, -1), round(end, -1), 10):
            if not count.get(chromosome + '_' + direction):
                count[chromosome + '_' + direction] = {}
            if not count[chromosome + '_' + direction].get(base):
                count[chromosome + '_' + direction][base] = number
            else:
                count[chromosome + '_' + direction][base] += number
    return count

def filter_isoforms(isoforms, count, isoform_names, chromosome, psl_info, overhangs):
    counter = 0
    for isoform in sorted(isoform_names):
        coverage_list = []
        info = psl_info[isoform]
        overhang5 = int(info[11])
        overhang3 = int(info[10]) - int(info[12])
        start = int(info[15])
        end = int(info[16])
        direction = info[8]
        number = int(isoform.split('_')[-1])
        for base in np.arange(round(start, -1), round(end, -1), 10):
            coverage_list.append(count[chromosome + '_' + direction][base])
        max_coverage = max(coverage_list)
        if number >= minimum_reads and (number/max_coverage) >= minimum_ratio:
             sequence = isoforms[isoform]
             if overhangs[0] <= overhang5 <= overhangs[1] and overhangs[2] <= overhang3 <= overhangs[3]:
                 out2.write('>%s\n%s\n' %(isoform, sequence))
                 out3.write('\t'.join(info) + '\n')

def look_for_contained_isoforms(psl_file, chromosome):
    filtered_isoforms = []
    covered = {}
    covered[chromosome] = {}
    psl_dict, psl_info = parse_clean_psl(psl_file, chromosome)
    for isoform, info in psl_dict.items():
        counter = 0
        coordinates = info[0]
        direction = info[1]
        start, end = coordinates[0], coordinates[-1]
        for terminus in [start, end]:
            for base in range(terminus-20, terminus+20, 1):
                if not covered[chromosome].get(base):
                    covered[chromosome][base] = set()
                covered[chromosome][base].add(isoform)

        for coord in coordinates:
            counter += 1
            if counter%2 == 1:
                start = coord
            if counter%2 == 0:
                end = coord
                for base in range(start-5, end+5, 1):
                    if not covered[chromosome].get(base):
                        covered[chromosome][base] = set()
                    covered[chromosome][base].add(isoform)

    for isoform, info in psl_dict.items():
        status = set()
        coordinates = info[0]
        direction = info[1]
        for isoform1 in psl_dict:
            status.add(isoform1)
        for coord in coordinates:
            counter += 1
            if counter%2 == 1:
                start = coord
            if counter%2 == 0:
                end = coord
                for base in range(start, end, 1):
                    status = status & covered[chromosome][base]

        minimum_overlap = len(status)

        if minimum_overlap == 1:
            filtered_isoforms.append(isoform)
        else:
            show = False
            isoform_abundance = int(isoform.split('_')[-1])
            if 'Isoform_53959' in isoform:
                show = True
            decision = []
            for match in status:
                if show:
                    print(match)

                if match != isoform:
                    match_abundance = int(match.split('_')[-1])
                    match_coordinates = psl_dict[match][0]
                    duplicate_dict = {}
                    for coord in match_coordinates[1:-1]:
                        for base in range(coord-5, coord+5, 1):

                            duplicate_dict[base] = 1
                    start, end = match_coordinates[0], match_coordinates[-1]
                    for terminus in [start, end]:
                        for base in range(terminus-20, terminus+20, 1):
                            duplicate_dict[base] = 1
                    matches = []
                    for coord in coordinates[1:-1]:
                        if coord in duplicate_dict:
                            matches.append(coord)
                    if len(matches) != len(coordinates[1:-1]):
                        if show:
                            print('check', len(matches), len(coordinates[1:-1]))
                        continue

                    matches = []
                    for coord in coordinates:
                        if coord in duplicate_dict:
                            matches.append(coord)
                    if show:
                        print(len(matches), len(coordinates))
                    if len(matches) == len(coordinates):
                        if isoform_abundance < match_abundance:
                            decision.append(False)
                            if show:
                                print('duplicate')

                    if (isoform_abundance/match_abundance) < internal_ratio:
                         if show:
                             print(match, isoform_abundance, match_abundance)
                         decision.append(False)

            if show:
                print(decision)
            if False not in decision:
                filtered_isoforms.append(isoform)

    return filtered_isoforms, psl_dict, psl_info

def reverse_complement(sequence):
    '''Returns the reverse complement of a sequence'''
    bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
    return ''.join([bases[x] for x in list(sequence)])[::-1]

def read_splice_file(SS_file):
    splice_dict = {}
    for line in open(SS_file):
        a = line.strip().split('\t')
        chromosome = a[0]
        start = int(a[1])
        end = int(a[2])
        name = a[3].split('_')[0]
        splice_dict[name] = (start, end)
    return splice_dict

def simplify(infile, outfile, namefile):
    isoforms = read_fasta(infile)
    out = open(outfile, 'w')
    out1 = open(namefile, 'w')
    counter = 0
    for isoform, sequence in isoforms.items():
        if len(sequence) > 0:
            counter += 1
            out.write('>Isoform' + '_'
                      + str(counter) + '_' + isoform.split('_')[-1]
                      + '\n' + sequence + '\n')
            out1.write(isoform + '\tIsoform' + '_' + str(counter) + '\n')
    out.close()
    out1.close()

def parse_clean_psl(psl_file, target_chromosome):
    psl_dict = {}
    psl_info = {}
    for line in open(psl_file):
        a = line.strip().split('\t')

        chromosome = a[13]
        if chromosome == target_chromosome:
            start = int(a[15])
            end = int(a[16])
            direction = a[8]
            name = a[9]
            psl_info[name] = a
            begin, span = int(a[15]), int(a[16])
            blocksizes = a[18].split(',')[:-1]
            blockstarts = a[20].split(',')[:-1]
            readstarts = a[19].split(',')[:-1]
            psl_dict[name] = [[],direction]
            for index in np.arange(0, len(blocksizes), 1):
                blockstart = int(blockstarts[index])
                blocksize = int(blocksizes[index])
                blockend = blockstart + blocksize
                psl_dict[name][0].append(blockstart)
                psl_dict[name][0].append(blockend)

    return psl_dict, psl_info

def collect_chromosomes(isoform_psl):
    chromosomes = set()
    for line in open(isoform_psl):
        a = line.strip().split('\t')
        chromosome = a[13]
        chromosomes.add(chromosome)
    chromosomes = sorted(list(chromosomes))
    return chromosomes

def main(infile):
    print('simplifying isoform names')
    temp_fasta = path + '/isoform_tmp.fasta'
    simplify(infile, temp_fasta, path + '/Isoform_long_names.txt')

    os.system('%s -i %s -a %s -o %s -c %s -e %s'% ('python3 postprocessingIsoforms.py', temp_fasta,
                                                   adapter, path, configIn, ends))
    print('reading fasta')
    processed_isoforms = path + 'Isoforms_full_length_consensus_reads.fasta'
    isoforms = read_fasta(processed_isoforms)
    sam_file = path + '/Isoforms.aligned.out.sam'
    psl_file = path + '/Isoforms.aligned.out.psl'
    clean_psl_file = path + '/Isoforms.aligned.out.clean.psl'
    print('aligning reads')
    os.system('%s --secondary=no -ax splice -t %s %s %s > %s ' %(minimap2, minimap2_threads, genome, processed_isoforms, sam_file))
    os.system('%s -i %s > %s ' %(emtrey, sam_file, psl_file))
    os.system('%s %s %s ' %('python3 clean_psl.py', psl_file, clean_psl_file))
    print('collecting chromosomes')
    chromosomes = collect_chromosomes(clean_psl_file)
    for chromosome in chromosomes:
        print('finding fully contained isoforms on', chromosome)
        isoform_names, psl_dict, psl_info = look_for_contained_isoforms(clean_psl_file, chromosome)
        print('getting isoform loci read counts')
        count = get_count(isoform_names, chromosome, psl_dict)
        print('final filtering step')
        filter_isoforms(isoforms, count, isoform_names, chromosome, psl_info, overhangs)

main(infile)

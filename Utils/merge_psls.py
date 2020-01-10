#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import os
import sys
import numpy as np

upstream_buffer = 10
downstream_buffer = 50
minimum_read_count = 1

inFile = sys.argv[1]

def find_peaks(starts, ends):
    start_peaks, end_peaks = {}, {}
    for position in sorted(starts):
        if list(starts).count(position) >= minimum_read_count:
            if not start_peaks.get(position):
                for shift in np.arange(-upstream_buffer, downstream_buffer):
                    start_peaks[position+shift] = position
    for position in sorted(ends, reverse=True):
        if list(ends).count(position) >= minimum_read_count:
            if not end_peaks.get(position):
                for shift in np.arange(-downstream_buffer, upstream_buffer, 1):
                    end_peaks[position+shift] = position
    return start_peaks, end_peaks

def get_alignment_direction(sam_file):
    direction_dict = {}
    direction_conversion = {}
    direction_conversion['+'] = '+'
    direction_conversion['-'] = '-'

    for line in open(sam_file):
        if line[0] != '@':
            a = line.strip().split('\t')
            read_name = a[0]
            for entry in a[10:]:
                if 'XS:A:' in entry:
                    direction = entry.split('XS:A:')[1]
                    direction_dict[read_name] = direction_conversion[direction]
    return direction_dict

def sort_reads_into_splice_junctions(infile, direction_dict, start_end_dict, cell, out_all, rep, redirect_dict):
    for line in open(infile):
        out_all.write(line)
        a = line.strip().split('\t')
        read_chromosome = a[13]
        name = a[9]
        if name in direction_dict:
            read_direction = direction_dict[name] # ignores read direction
        else:
            read_direction = '+'
        start, end = int(a[15]), int(a[16])

        identity = read_chromosome + '_' + read_direction + '_'

        blocksizes = a[18].split(',')[:-1]
        blockstarts = a[20].split(',')[:-1]
        readstarts = a[19].split(',')[:-1]

        for x in range(0, len(blocksizes)-1):
            blockstart = int(blockstarts[x])
            blocksize = int(blocksizes[x])
            left_splice_site = blockstart + blocksize
            right_splice_site = int(blockstarts[x+1])
            if 61013896 < right_splice_site < 61013941:
                print('before', right_splice_site)
            if left_splice_site in redirect_dict[read_chromosome][read_direction]['left']:
                left_splice_site = redirect_dict[read_chromosome][read_direction]['left'][left_splice_site]
            if right_splice_site in redirect_dict[read_chromosome][read_direction]['right']:
                right_splice_site = redirect_dict[read_chromosome][read_direction]['right'][right_splice_site]
            if 61013896 < right_splice_site < 61013941:
                print('after', right_splice_site)

            identity += str(left_splice_site) + '-' \
                            + str(right_splice_site) + '~'

        if not start_end_dict.get(identity):
            start_end_dict[identity] = set()
        start_end_dict[identity].add((start, end, tuple(a), read_direction, cell, rep))
    return start_end_dict

def combine_equivalent_isoforms(positions, count_dict, pos_dict, interval):
    merged_positions = set()
    for position in sorted(positions, key=lambda x:(round(x[0], -1), -round(x[1], -1))):

        start, end, a, read_direction, cell, rep = position[0], position[1], position[2], position[3], position[4], position[5]
        matched = False
        for x in range(start-interval, start+interval):
            if matched:
                pass
            else:
                for y in range(end-interval, end+interval):
                    if x in pos_dict and y in pos_dict[x]:
                        matched_position = pos_dict[x][y]
                        matched = True
        if matched:
            if position != matched_position:
                count_dict[matched_position[2][9]].add(rep)
            merged_positions.add(matched_position)
            for x in range(start-interval, start+interval):
                if x not in pos_dict:
                    pos_dict[x] = {}
                for y in range(end-interval, end+interval):
                    pos_dict[x][y] = matched_position
        else:
            count_dict[a[9]] = set()
            count_dict[a[9]].add(rep)
            merged_positions.add(position)
            for x in range(start-interval, start+interval):
                if x not in pos_dict:
                    pos_dict[x] = {}
                for y in range(end-interval, end+interval):
                    pos_dict[x][y] = position
    return merged_positions, count_dict, pos_dict

def define_start_end_sites(start_end_dict, out, read_dict, out_fasta):
    isoform_counter, isoform_dict, subread_pointer = 0, {}, {}
    number_of_isoforms = len(start_end_dict)
    counter = 0
    count_dict = {}
    interval = 10
    for identity in start_end_dict:
        positions = start_end_dict[identity]
        pos_dict = {}

        pre = len(positions)
        positions, count_dict, pos_dict = combine_equivalent_isoforms(positions, count_dict, pos_dict, interval)
        positions, count_dict, pos_dict = combine_equivalent_isoforms(positions, count_dict, pos_dict, interval)
        positions, count_dict, pos_dict = combine_equivalent_isoforms(positions, count_dict, pos_dict, interval)
        post = len(positions)

        for position in positions:
            a = list(position[2])
            cell = position[4]
            rep = position[5]
            abundance = len(count_dict[a[9]])
            write = True
            if 'Cluster' in rep and abundance>0:
                write = True
            if abundance >= 1:
                write = True
            if write:
                a[9] = rep + '-' + cell + '-' + a[9] + '-' + str(abundance)
                out.write(('\t').join(a) + '\n')

def read_fasta(inFile, rep, cell, readDict):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    lastHead = ''
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            if lastHead:
                readDict[lastHead] = ''.join(readDict[lastHead])
            lastHead = rep + '-' + cell + '-' + ('_').join(line[1:].split('_')[:-1])
            readDict[lastHead] = []
        else:
            readDict[lastHead].append(line.upper())
    if readDict:
        readDict[lastHead] = ''.join(readDict[lastHead])
    return readDict

def define_splice_junctions(infile, direction_dict, cell, rep, splice_junction_dict):
    for line in open(infile):
        a = line.strip().split('\t')
        read_chromosome = a[13]
        if read_chromosome not in splice_junction_dict:
            splice_junction_dict[read_chromosome] = {}
            splice_junction_dict[read_chromosome]['+'] = {}
            splice_junction_dict[read_chromosome]['-'] = {}
            splice_junction_dict[read_chromosome]['+']['left'] = {}
            splice_junction_dict[read_chromosome]['-']['left'] = {}
            splice_junction_dict[read_chromosome]['+']['right'] = {}
            splice_junction_dict[read_chromosome]['-']['right'] = {}

        name = a[9]
        if name in direction_dict:
            read_direction = direction_dict[name] # ignores read direction
        else:
            read_direction = a[8]#'+'
        start, end = int(a[15]), int(a[16])

        blocksizes = a[18].split(',')[:-1]
        blockstarts = a[20].split(',')[:-1]
        readstarts = a[19].split(',')[:-1]

        for x in range(0, len(blocksizes)-1):
            blockstart = int(blockstarts[x])
            blocksize = int(blocksizes[x])
            left_splice_site = blockstart + blocksize
            right_splice_site = int(blockstarts[x+1])
            if not left_splice_site in splice_junction_dict[read_chromosome][read_direction]['left']:
                splice_junction_dict[read_chromosome][read_direction]['left'][left_splice_site] = 0
            splice_junction_dict[read_chromosome][read_direction]['left'][left_splice_site] += 1
            if not right_splice_site in splice_junction_dict[read_chromosome][read_direction]['right']:
                splice_junction_dict[read_chromosome][read_direction]['right'][right_splice_site] = 0
            splice_junction_dict[read_chromosome][read_direction]['right'][right_splice_site] += 1
    return splice_junction_dict

def populate_redirect_dict(splice_junction_dict):
    redirect_dict = {}
    for chromosome in splice_junction_dict:
        redirect_dict[chromosome] = {}
        for direction in ['+', '-']:
            redirect_dict[chromosome][direction] = {}
            for orientation in ['left', 'right']:
                redirect_dict[chromosome][direction][orientation] = {}
                for site, abundance in splice_junction_dict[chromosome][direction][orientation].items():
                    if abundance > 20:
                        for replacement in (site-1, site+1):
                            if replacement in splice_junction_dict[chromosome][direction][orientation]:
                                previous = splice_junction_dict[chromosome][direction][orientation][replacement]
                                if previous<abundance/10:
                                    redirect_dict[chromosome][direction][orientation][replacement] = site
                            else:
                                redirect_dict[chromosome][direction][orientation][replacement] = site
    return redirect_dict

def main():
    print('defining spliceforms')
    out = open('Mandalorion_FLAIR_merged.psl', 'w')
    out_fasta = open('Mandalorion_FLAIR_merged.fasta', 'w')
    out_all = open('Mandalorion_FLAIR_all.psl', 'w')
    start_end_dict = {}
    read_dict = {}
    splice_junction_dict = {}

    for line in open(inFile):
        a = line.strip().split('\t')
        psl = a[0]
        fasta = a[1]
        cell = a[2]
        rep = a[3]
        direction_dict1 = {}
        splice_junction_dict = define_splice_junctions(psl, direction_dict1, cell, rep, splice_junction_dict)

    redirect_dict = populate_redirect_dict(splice_junction_dict)

    for line in open(inFile):
        a = line.strip().split('\t')
        psl = a[0]
        fasta = a[1]
        cell = a[2]
        rep = a[3]
        direction_dict1 = {}
        start_end_dict = sort_reads_into_splice_junctions(psl, direction_dict1, start_end_dict, cell, out_all, rep, redirect_dict)
        read_dict = read_fasta(fasta, rep, cell, read_dict)

    print('merging equivalent transcripts')
    define_start_end_sites(start_end_dict, out, read_dict, out_fasta)

if __name__ == '__main__':
    main()

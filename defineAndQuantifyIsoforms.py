#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import os
import sys
import numpy as np

path = sys.argv[2]
infile = sys.argv[1]
upstream_buffer = int(sys.argv[4])
downstream_buffer = int(sys.argv[3])
subreads=sys.argv[5]
fasta_file=sys.argv[6]

minimum_read_count = 3

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

def collect_splice_events(path):
    splice_dict = {}
    for line in open(path + '/SS.bed'):
        a = line.strip().split('\t')
        chromosome = a[0]
        if not splice_dict.get(chromosome):
            splice_dict[chromosome] = {}
        splice_left =int(a[1])
        splice_right = int(a[2])
        for base in np.arange(splice_left, splice_right+1):
             splice_dict[chromosome][base] = a[3].split('_')[0]
    return splice_dict

def sort_reads_into_splice_junctions(splice_dict,
                                     fasta_file, infile):
    start_end_dict, readDict = {}, {}
    tempSeqs, headers, sequences = [], [], []
    for line in open(fasta_file):
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
        readDict[headers[i].split('_')[0]] = [headers[i], sequences[i]]
    read_dict = readDict

    for line in open(infile):
        a = line.strip().split('\t')
        read_chromosome, read_direction = a[13], a[8]
        name = a[9].split('_')[0]
        read_direction = '+' # ignores read direction
        start, end = int(a[15]), int(a[16])
        if end - start > 500:
            if read_direction == '+':
                left_extra, right_extra = int(a[11]), int(a[10]) - int(a[12])
            if read_direction == '-':
                right_extra, left_extra = int(a[11]), int(a[10]) - int(a[12])

            failed = False
            identity = read_chromosome + '_'

            blocksizes = a[18].split(',')[:-1]
            blockstarts = a[20].split(',')[:-1]
            readstarts = a[19].split(',')[:-1]

            for x in range(0, len(blocksizes)-1):
                blockstart = int(blockstarts[x])
                blocksize = int(blocksizes[x])
                left_splice = blockstart + blocksize
                right_splice = int(blockstarts[x+1])
                if right_splice - left_splice > 50:
                    if not splice_dict.get(read_chromosome):
                        failed = True
                        break
                    else:
                        if not splice_dict[read_chromosome].get(left_splice) \
                                or not splice_dict[read_chromosome].get(right_splice):
                            failed = True
                            break
                    left_splice_site = splice_dict[read_chromosome][left_splice]
                    right_splice_site = splice_dict[read_chromosome][right_splice]
                    identity += str(left_splice_site) + '-' \
                                + str(right_splice_site) + '~'
            if not failed:
                if not start_end_dict.get(identity):
                    start_end_dict[identity] = []
                start_end_dict[identity].append((start, end,
                                                 '>' + read_dict[name][0] + '\n'
                                                 + read_dict[name][1] + '\n',
                                                 left_extra,
                                                 right_extra,
                                                 read_direction))
    return start_end_dict

def define_start_end_sites(start_end_dict, individual_path):
    left_extras, right_extras = {}, {}
    file_set = set()
    isoform_counter, isoform_dict,subread_pointer = 0, {}, {}
    number_of_isoforms = len(start_end_dict)
    counter = 0
    for identity in start_end_dict:
        counter += 1
        print('processing spliceform ', counter, 'of ',
              number_of_isoforms, ' into isoforms')
        starts = []
        ends = []
        positions = start_end_dict[identity]
        length_of_positions = len(positions)
        indexes = np.random.choice(range(0, length_of_positions, 1),
                                   size=min(10000, length_of_positions),
                                   replace=False)
        new_positions = []
        for index in indexes:
            new_positions.append(positions[index])
        positions = new_positions
        for position in positions:
            starts.append(int(position[0]))
            ends.append(int(position[1]))

        start_dict, end_dict = find_peaks(starts, ends)
        matched_positions = []
        combination_counts = {}
        left_extras[identity], right_extras[identity] = {}, {}

        for start, end, read, left_extra, right_extra, read_direction in positions:
            try:
                left = start_dict[int(start)]
                right = end_dict[int(end)]
                if not left_extras[identity].get((left, right)):
                    left_extras[identity][(left,right)] = []
                    right_extras[identity][(left,right)] = []

                left_extras[identity][(left, right)].append(int(left_extra))
                right_extras[identity][(left, right)].append(int(right_extra))

                if not combination_counts.get((left, right)):
                    combination_counts[(left, right)] = 1
                else:
                    combination_counts[(left, right)] += 1
                matched_positions.append((left, right, read, read_direction))
            except:
                pass
        for left, right, read, read_direction in matched_positions:
            medianLeft = np.median(left_extras[identity][(left, right)])
            medianRight = np.median(right_extras[identity][(left, right)])
            new_identity = identity + '_' + read_direction + '_' + str(left) \
                           + '_' + str(right) + '_' \
                           + str(round(medianLeft, 2)) \
                           + '_' + str(round(medianRight, 2))
            if not isoform_dict.get(new_identity):
                isoform_counter += 1
                isoform_dict[new_identity] = isoform_counter

            subfolder = str(int(isoform_dict[new_identity]/4000))
            if subfolder not in os.listdir(individual_path + '/parsed_reads/'):
                os.makedirs(individual_path + '/parsed_reads/' + subfolder)
            filename = subfolder + '/Isoform' + str(isoform_dict[new_identity])
            out_reads_fasta = open(individual_path + '/parsed_reads/'
                                   + filename + '.fasta', 'a')
            out_reads_subreads = open(individual_path + '/parsed_reads/'
                                      + filename + '_subreads.fastq', 'w')
            out_reads_fasta.write(read)
            out_reads_fasta.close()
            out_reads_subreads.close()

            file_set.add(individual_path + '/parsed_reads/' + filename
                         + '.fasta' + '\t' + individual_path
                         + '/parsed_reads/' + filename + '_subreads.fastq'
                         + '\t' + new_identity + '\n')
            read = read.split()[0].split('_')[0][1:]

            subread_pointer[read] = individual_path + '/parsed_reads/' \
                                    + filename + '_subreads.fastq'
            # for subread, sequence, qual in subread_list:
            #     out_reads_subreads.write(subread + '\n' + sequence
            #                              + '\n+\n' + qual + '\n')

    out = open(individual_path + 'isoform_list', 'a')
    for item in file_set:
        out.write(item)
    out.close()
    return subread_pointer

def read_subreads(seq_file, infile, subread_pointer):
    seq_file_open = open(seq_file, 'r')

    lineNum, lastPlus = 0, False
    for line in open(seq_file):
        line = line.rstrip()
        if not line:
            continue
        # make an entry as a list and append the header to that list
        if lineNum % 4 == 0 and line[0] == '@':
            if lastPlus:  # chrom_reads needs to contain root_names
                 filepath = subread_pointer.get(root_name)
                 if filepath:
                     outfile = open(filepath,'a')
                     outfile.write('@%s\n%s\n+\n%s\n' %(name, seq, qual))
                     outfile.close()
            name = line[1:]
            root_name = name.split('_')[0]

        if lineNum % 4 == 1:
            seq = line
        if lineNum % 4 == 2:
            lastPlus = True
        if lineNum % 4 == 3 and lastPlus:
            qual = line
        lineNum += 1

def main():
    individual_path = path
    os.system('mkdir ' + individual_path + '/parsed_reads')
    os.system('rm -r ' + individual_path + '/parsed_reads/*')
    out = open(individual_path + 'isoform_list', 'w')
    out.close()

    print('reading splice sites')
    splice_dict = splice_dict = collect_splice_events(path)
    print('defining spliceforms')
    start_end_dict = sort_reads_into_splice_junctions(splice_dict,
                                                      fasta_file, infile)
    print('defining isoforms')
    subread_pointer = define_start_end_sites(start_end_dict, individual_path)
    print('reading subreads')
    read_subreads(subreads, infile, subread_pointer)

if __name__ == '__main__':
    main()

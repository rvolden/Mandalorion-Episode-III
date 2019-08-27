#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import numpy as np

infile =  sys.argv[1]
out_path = sys.argv[2]
cutoff = float(sys.argv[3])
genome_file = sys.argv[4]
refine = sys.argv[5]
sam_file =  sys.argv[6]

minimum_read_count = 2

splice_site_width = 5

def scan_for_best_bin(entry, distance_range, iterator_shift, density_dict,
                      base_cutoff_min, base_cutoff_max,
                      peak_areas, chromosome, side):
    '''
    Stuff about the inputs, outputs, and how it works
    '''

    best_extra_list, peak_center, bases = [], 0, []
    coverage_area, best_direction_l, best_direction_r = [], [], []
    for x in distance_range:
        extra_list_bases, extra_list_expression = [], []
        direction_l, direction_r = {}, {}
        coverage_set = []
        called = False
        for y in distance_range:
            try:
               called = peak_areas[chromosome][side][entry+x+y]
            except KeyError:
               pass
        if not called:
             highest_y = 0
             highest_y_pos = 0
             for y in distance_range:
                try:
                    for item in density_dict[entry+x+y]:
                        extra_list_bases.append(item[0])
                        extra_list_expression.append(1)
                        if not direction_l.get(item[4]):
                            direction_l[item[4]] = 1
                        else:
                            direction_l[item[4]] += 1
                        if not direction_r.get(item[5]):
                            direction_r[item[5]] = 1
                        else:
                            direction_r[item[5]] += 1
                        for covered_position in item[3]:
                            coverage_set.append(covered_position)
                except:
                    pass

        if base_cutoff_min <= np.median(extra_list_bases) <= base_cutoff_max:
            if sum(extra_list_expression) > sum(best_extra_list):
                best_extra_list = extra_list_expression
                peak_center = entry + x
                bases = extra_list_bases
                coverage_area = coverage_set
                best_direction_l = direction_l
                best_direction_r = direction_r

    return best_extra_list, peak_center, bases, coverage_area, \
           best_direction_l, best_direction_r

def determine_coverage(coverage_area, chromosome, reverse,
                       peak_center, histo_coverage):
    coverage = [0]
    coverage_area2 = []
    for covered_position in set(coverage_area):
        if coverage_area.count(covered_position) > 1:
            coverage_area2.append(covered_position)

    coverage_area = sorted(coverage_area2, reverse=reverse)
    counter = 0
    for base_f in coverage_area:
        count = 0
        if not reverse:
            if base_f > peak_center:
                count = 1
        elif reverse:
            if base_f < peak_center:
                count = 1
        if count == 1:
            if counter <= 3:
                counter += 1
                base_f = myround(base_f)
                try:
                    coverage.append(histo_coverage[chromosome][base_f])
                except KeyError:
                    pass
            else:
                break
    coverage = max(coverage)
    return coverage, coverage_area

def myround(x, base=10):
    '''Rounds to the nearest base'''
    return int(base * round(float(x)/base))

def find_peaks(density_dict, out, peaks, reverse, cutoff, base_cutoff_min,
               base_cutoff_max, histo_coverage, side, peak_areas,chromosome):
    '''
    Insert docstring
    '''

    if not reverse:
        distance_range = range(-splice_site_width, splice_site_width)
        iterator_shift = 1
    if reverse:
        distance_range = range(splice_site_width, -splice_site_width, -1)
        iterator_shift =- 1

    entry_list = []
    for entry in density_dict:
        entry_list.append([entry, density_dict[entry]])
    for entry, density in sorted(entry_list,
                                 key=lambda x: len(x[1]),
                                 reverse=True):

        if len(density) >= minimum_read_count:
            if not peak_areas[chromosome][side].get(entry):

                best_extra_list, peak_center, bases, \
                coverage_area, best_direction_l, best_direction_r \
                = scan_for_best_bin(entry, distance_range, iterator_shift,
                                    density_dict, base_cutoff_min,
                                    base_cutoff_max, peak_areas,
                                    chromosome, side)

                coverage, coverage_area \
                = determine_coverage(coverage_area, chromosome, reverse,
                                     peak_center, histo_coverage)


                if coverage > 0:
                    proportion = round(sum(best_extra_list)/coverage, 3)
                    if proportion > cutoff:
                        if 'TSS' in best_direction_l:
                            Left_TSS = best_direction_l['TSS']
                        else:
                            Left_TSS = 0
                        if 'TES' in best_direction_l:
                            Left_TES = best_direction_l['TES']
                        else:
                            Left_TES = 0
                        if 'TSS' in best_direction_r:
                            Right_TSS = best_direction_r['TSS']
                        else:
                            Right_TSS = 0
                        if 'TES' in best_direction_r:
                            Right_TES = best_direction_r['TES']
                        else:
                            Right_TES = 0
                        Left_to_Right = Left_TSS + Right_TES
                        Right_to_Left = Left_TES + Right_TSS
                        Type = '-'

                        if Left_to_Right < Right_to_Left and reverse:
                            Type = '3'
                        elif Left_to_Right < Right_to_Left and not reverse:
                            Type = '5'
                        elif Left_to_Right > Right_to_Left and reverse:
                            Type = '5'
                        elif Left_to_Right > Right_to_Left and not reverse:
                            Type = '3'

                        if Type != '-':
                            peaks += 1
                            out.write(chromosome + '\t'
                                      + str(peak_center-splice_site_width)
                                      + '\t' + str(peak_center+splice_site_width)
                                      + '\t' + str(Type) + side + str(peaks) + '_'
                                      + str(peak_center-splice_site_width) + '_'
                                      + str(peak_center+splice_site_width) + '_'
                                      + str(proportion) + '\t' + str(peaks) + '\n')
                            for base in range(peak_center - splice_site_width,
                                              peak_center + splice_site_width):
                                peak_areas[chromosome][side][base] = 1
        else:
            break

    return peaks, peak_areas

def collect_reads(infile, sam_file, chromosome_list):
    histo_left_bases, histo_right_bases = {}, {}
    chromosome_list_left, chromosome_list_right = chromosome_list, chromosome_list
    for chromosome in chromosome_list:
        histo_left_bases[chromosome] = {}
        histo_right_bases[chromosome] = {}


    histo_coverage = {}
    base_cutoff_min, base_cutoff_max = 0, 5

    total = 0
    length = 0

    direction_dict = get_alignment_direction(sam_file)
    for line in open(infile):
            total += 1
            a = line.strip().split('\t')
            chromosome = a[13]
            name = a[9].split('_')[0]
            direction = a[8]
            if direction_dict.get(name):
                direction = direction_dict[name]
            if not histo_coverage.get(chromosome):
                histo_coverage[chromosome] = {}
            score, name = int(a[0]), a[9]
            length = int(a[10])
            begin, span = int(a[15]), int(a[16])
            blocksizes = a[18].split(',')[:-1]
            blockstarts = a[20].split(',')[:-1]
            readstarts = a[19].split(',')[:-1]

            if direction == '+':
                start_seq, end_seq = 'S', 'E'
                left_match, right_match = 'TSS', 'TES'

            else:
                start_seq, end_seq= 'E', 'S'
                left_match, right_match = 'TES', 'TSS'

            coverage_set = set()
            low_bounds, up_bounds = [], []
            aligned_bases = 0
            for x in range(0, len(blocksizes)):
                blockstart = int(blockstarts[x])
                blocksize = int(blocksizes[x])
                readstart = int(readstarts[x])
                aligned_bases += blocksize
                blockend = blockstart + blocksize
                for y in range(0, blocksize, 10):
                    rounded = myround(blockstart + y)
                    coverage_set.add(rounded)
                for yy in range(y, blocksize):
                    rounded = myround(blockstart + yy)
                    coverage_set.add(rounded)
                if blockstart != begin:
                    up_bounds.append(blockstart)
                if blockend != span:
                    low_bounds.append(blockend)

            for rounded in coverage_set:
                try:
                    histo_coverage[chromosome][rounded] += 1
                except:
                    histo_coverage[chromosome][rounded] = 1

            if aligned_bases/length > 0.70:
                for low_bound in low_bounds:
                    chromosome_list_left.add(chromosome)
                    if not histo_left_bases.get(chromosome):
                        histo_left_bases[chromosome] = {}
                    if not histo_left_bases[chromosome].get(low_bound):
                        histo_left_bases[chromosome][low_bound] = []
                    histo_left_bases[chromosome][low_bound].append([0, begin,
                                                                    span, coverage_set,
                                                                    left_match,
                                                                    right_match])
                for up_bound in up_bounds:
                    chromosome_list_right.add(chromosome)
                    if not histo_right_bases.get(chromosome):
                        histo_right_bases[chromosome] = {}
                    if not histo_right_bases[chromosome].get(up_bound):
                        histo_right_bases[chromosome][up_bound] = []
                    histo_right_bases[chromosome][up_bound].append([0, begin,
                                                                    span, coverage_set,
                                                                    left_match,
                                                                    right_match])

    chromosome_list = chromosome_list_left & chromosome_list_right
    return histo_left_bases, histo_right_bases, chromosome_list, histo_coverage

def parse_genome(input_file, left_bounds, right_bounds):
    chromosome_list = set()
    gene_dict = {}
    for line in open(input_file):
        a = line.strip().split('\t')
        if len(a) > 7:
             if a[2] == 'exon':
                 testKey = a[8].split('transcript_id "')[1].split('"')[0]
                 if not gene_dict.get(testKey):
                     gene_dict[testKey] = []
                 gene_dict[testKey].append((a[0], a[3], a[4], a[6]))

    read_list = []
    for transcript_id in gene_dict:
        transcript_data = gene_dict[transcript_id]

        chromosome = transcript_data[0][0]
        chromosome_list.add(chromosome)
        if not right_bounds.get(chromosome):
            left_bounds[chromosome], right_bounds[chromosome]= {}, {}
            left_bounds[chromosome]['5'], right_bounds[chromosome]['5'] = [], []
            left_bounds[chromosome]['3'], right_bounds[chromosome]['3'] = [], []

        start = sorted(transcript_data, key=lambda x: int(x[1]))[0][1]
        end = sorted(transcript_data, key=lambda x: int(x[2]), reverse=True)[0][2]

        for entry in transcript_data:
            if entry[1] != start:
                if entry[3] == '+':
                    right_bounds[chromosome]['3'].append(int(entry[1])-1)
                elif entry[3] == '-':
                    right_bounds[chromosome]['5'].append(int(entry[1])-1)
            if entry[2] != end:
                if entry[3] == '+':
                    left_bounds[chromosome]['5'].append(int(entry[2]))
                if entry[3] == '-':
                    left_bounds[chromosome]['3'].append(int(entry[2]))
    return chromosome_list, left_bounds, right_bounds

def make_genome_bins(bounds, side, peaks, chromosome, peak_areas):
    for type1 in ['5', '3']:
        covered = {}
        position_list = sorted(bounds[type1], key=int)
        for index1 in range(0, len(position_list)):
            if not covered.get(index1):
                sub_list = []
                sub_list.append(position_list[index1])
                for index2 in range(index1, len(position_list)):
                    if position_list[index2] - max(sub_list) <= splice_site_width:
                        sub_list.append(position_list[index2])
                        covered[index2] = 1
                    else:
                        break
                single = 0
                if len(sub_list) > 1:
                    splice_distances = []
                    for splice_pos in range(0, len(sub_list)-1):
                        splice_distances.append(int(sub_list[splice_pos+1])
                                                - int(sub_list[splice_pos]))
                    if min(splice_distances) > 3:
                        for x in range(0, len(sub_list), 1):
                            if x != 0:
                                start = int(sub_list[x]
                                            - ((sub_list[x] - sub_list[x-1])/2))
                            else:
                                start = int(sub_list[x]) - 1
                            if x != len(sub_list) - 1:
                                end = int(sub_list[x]
                                          + ((sub_list[x+1] - sub_list[x])/2))
                            else:
                                end = int(sub_list[x]) + 1

                            out.write(chromosome + '\t' + str(start) + '\t'
                                      + str(end) + '\t' + type1 + side
                                      + str(peaks) + '_' + str(start) + '_'
                                      + str(end) + '_A' + '\t' + str(peaks) + '\n')
                            for base in range(start, end):
                                peak_areas[chromosome][side][base] = 1
                            peaks += 1
                    else:
                         single = 1
                else:
                    single = 1
                if single == 1:
                    start = min(sub_list) - splice_site_width
                    end = max(sub_list) + splice_site_width
                    out.write(chromosome + '\t' + str(start) + '\t' + str(end)
                              + '\t' + type1 + side + str(peaks) + '_'
                              + str(start) + '_' + str(end) + '_A'
                              + '\t' + str(peaks) + '\n')
                    for base in range(start-1, end+1):
                        peak_areas[chromosome][side][base] = 1
                    peaks += 1

    return peaks, peak_areas

def get_alignment_direction(sam_file):
    direction_dict = {}
    direction_conversion = {}
    direction_conversion[('0','+')] = '+'
    direction_conversion[('0','-')] = '-'
    direction_conversion[('16','+')] = '-'
    direction_conversion[('16','-')] = '+'

    for line in open(sam_file):
        if line[0] != '@':
            a = line.strip().split('\t')
            read_name = a[0].split('_')[0]
            read_direction = a[1]
            if read_direction in ['0', '16']:
                for entry in a[10:]:
                    if 'ts:A:' in entry:
                        direction = entry.split('ts:A:')[1]
                        direction_dict[read_name] = direction_conversion[(read_direction, direction)]
    return direction_dict

left_bounds = {}
right_bounds = {}
print('parsing annotated splice sites')
chromosome_list,left_bounds, right_bounds = parse_genome(genome_file, left_bounds, right_bounds)
print(sorted(list(chromosome_list)))
Left_Peaks = 0
Right_Peaks = 0
print('collecting reads')
histo_left_bases, histo_right_bases, \
chromosome_list, histo_coverage = collect_reads(infile,sam_file,chromosome_list)
out = open(out_path + '/SS.bed', 'w')

peak_areas = {}
print(sorted(list(chromosome_list)))
for chromosome in sorted(list(chromosome_list)):
    peak_areas[chromosome] = {}
    peak_areas[chromosome]['l'] = {}
    peak_areas[chromosome]['r'] = {}
    if not left_bounds.get(chromosome):
        left_bounds[chromosome]={}
        left_bounds[chromosome]['5'], left_bounds[chromosome]['3'] = [], []
    if not right_bounds.get(chromosome):
        right_bounds[chromosome] = {}
        right_bounds[chromosome]['5'], right_bounds[chromosome]['3'] = [], []
    print('now processing ', chromosome)
    if 'g' in refine:
        Left_Peaks_old = Left_Peaks
        Right_Peaks_old = Right_Peaks
        Left_Peaks, peak_areas = make_genome_bins(left_bounds[chromosome], 'l',
                                                  Left_Peaks, chromosome, peak_areas)
        Right_Peaks, peak_areas = make_genome_bins(right_bounds[chromosome], 'r',
                                                   Right_Peaks, chromosome, peak_areas)
        print('Parsed annotation-based splice-sites',
              Left_Peaks - Left_Peaks_old,
              Right_Peaks - Right_Peaks_old)
    Left_Peaks_old = Left_Peaks
    Right_Peaks_old = Right_Peaks
    Left_Peaks, peak_areas = find_peaks(histo_left_bases[chromosome],
                                        out, Left_Peaks, True, cutoff,
                                        0, 5, histo_coverage, 'l',
                                        peak_areas, chromosome)
    Right_Peaks, peak_areas = find_peaks(histo_right_bases[chromosome],
                                         out, Right_Peaks, False, cutoff,
                                         0, 5, histo_coverage, 'r',
                                         peak_areas, chromosome)
    print('Detected read-based splice-sites',
          Left_Peaks - Left_Peaks_old,
          Right_Peaks - Right_Peaks_old)

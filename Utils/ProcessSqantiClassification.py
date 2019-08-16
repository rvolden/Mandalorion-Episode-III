import sys
import os
import math
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--sqanti_classification_file', type=str)
parser.add_argument('-p', '--isoform_psl_file', type=str)
parser.add_argument('-f', '--isoform_fasta_file', type=str)
parser.add_argument('-o', '--output_folder', type=str)
parser.add_argument('-g', '--gtf_file', type=str)

args = parser.parse_args()
classification = args.sqanti_classification_file
psl_file = args.isoform_psl_file
fasta_file = args.isoform_fasta_file
path = args.output_folder
gtf_file = args.gtf_file

def read_genome_gtf(gtf_file, target_chromosome):
    gene_dict, exon_dict, terminal_dict, transcript_dict = {}, {}, {}, {}
    for line in open(gtf_file):
        a = line.strip().split('\t')
        if len(a) > 7 and a[0] == target_chromosome:
            if a[2] == 'exon':
                testKey = a[8].split('transcript_id "')[1].split('"')[0]
                if not transcript_dict.get(testKey):
                    transcript_dict[testKey] = []
                transcript_dict[testKey].append((a[0], a[3], a[4], a[6]))

    read_list = []
    for transcript_id in transcript_dict:
        transcript_data = transcript_dict[transcript_id]
        chromosome = transcript_data[0][0]
        if not gene_dict.get(chromosome):
            gene_dict[chromosome] = {}
            exon_dict[chromosome] = {}

        start = sorted(transcript_data, key=lambda x: int(x[1]))[0][1]
        end = sorted(transcript_data, key=lambda x: int(x[2]), reverse=True)[0][2]

        for entry in range(int(start), int(end), 1):
            gene_dict[chromosome][entry] = transcript_id

        for exon in transcript_data:
            for entry in range(int(exon[1]), int(exon[2]), 1):
                exon_dict[chromosome][entry] = transcript_id

        if not terminal_dict.get(chromosome):
            terminal_dict[chromosome] = {}
        for x in range(int(start)-10, int(start)+10, 1):
            terminal_dict[chromosome][x] = transcript_id
        for x in range(int(end)-10, int(end)+10, 1):
            terminal_dict[chromosome][x] = transcript_id

    return gene_dict, exon_dict, terminal_dict

def read_contigs(contig_file, gene_dict, exon_dict, target_chromosome):
    terminal_dict = {}
    for line in open(contig_file):
        a = line.strip().split('\t')
        chromosome = a[13]
        if chromosome == target_chromosome:
            start = int(a[15])
            end = int(a[16])

            if not terminal_dict.get(chromosome):
                terminal_dict[chromosome] = {}
            for x in range(start-10, start+10, 1):
                terminal_dict[chromosome][x] = a[9]
            for x in range(end-10, end+10, 1):
                terminal_dict[chromosome][x] = a[9]

            if not gene_dict.get(chromosome):
                gene_dict[chromosome] = {}
                exon_dict[chromosome] = {}

            begin, span = int(a[15]), int(a[16])
            blocksizes = a[18].split(',')[:-1]
            blockstarts = a[20].split(',')[:-1]
            readstarts = a[19].split(',')[:-1]
            for entry in range(begin, span, 1):
                if not gene_dict[chromosome].get(entry):
                    gene_dict[chromosome][entry] = a[9]

            for x in range(0, len(blocksizes), 1):
                blockstart = int(blockstarts[x])
                blocksize = int(blocksizes[x])
                blockend = blockstart + blocksize
                for entry in range(blockstart, blockend):
                    exon_dict[chromosome][entry] = a[9]

    return gene_dict,exon_dict, terminal_dict

def test_overlap(name, exonstart, exonend, exon_dict, direction, chromosome, new_exons):
    exon_overlap = False
    background_set = set([1])
    overlap_set = set()
    dec_set = set()
    for entry in range(exonstart, exonend, 1):
        if exon_dict[chromosome].get(entry):
            exon_overlap = True
            overlap_set.add(exon_dict[chromosome][entry])
        else:
            exon_overlap = False
        dec_set.add(exon_overlap)

    if False in dec_set and not overlap_set - background_set:
        new_exons.write(chromosome + '\t' + str(exonstart) + '\t' \
                        + str(exonend) + '\t' + name + '\t0\t' + direction + '\n')
        for entry in range(exonstart, exonend, 1):
            exon_dict[chromosome][entry] = 1
    return exon_dict

def read_isoforms(
        loci_content, gene_dict, exon_dict, terminal_dict,
        target_chromosome, targets, all_loci, new_genes,
        new_exons, new_TSS, new_polyA, all_loci_verbose, count_loci):

    for locus in loci_content:
        starts, ends, directions = [], [], []
        chromosome = locus.split('_')[0]
        for isoform in loci_content[locus]:
            starts.append(isoform[0])
            ends.append(isoform[1])
            directions.append(isoform[2])

        locus_start = min(starts)
        locus_end = max(ends)
        if chromosome == target_chromosome:
            if not terminal_dict.get(chromosome):
                terminal_dict[chromosome] = {}
            if not gene_dict.get(chromosome):
                gene_dict[chromosome] = {}
                exon_dict[chromosome] = {}

            all_loci.write(locus + '\t' + str(locus_start) + '\t' \
                           + str(locus_end) + '\t' \
                           + str(len(loci_content[locus])) + '\n')
            all_loci_verbose.write(locus + '\t' + str(locus_start) \
                                   + '\t' + str(locus_end) + '\t')
            for isoform in loci_content[locus]:
                all_loci_verbose.write(isoform[6] + ',')
            all_loci_verbose.write('\n')

            gene_match = {}
            gene_match['-'] = 1
            gene_overlap = False
            for entry in range(locus_start, locus_end, 1):
                    if gene_dict[chromosome].get(entry):
                        gene_overlap = True
                        if not gene_match.get(gene_dict[chromosome][entry]):
                            gene_match[gene_dict[chromosome][entry]] = 0
                        gene_match[gene_dict[chromosome][entry]] += 1

            match_list = []
            for match in gene_match:
                match_list.append((match, gene_match[match]))
            sorted_list = sorted(match_list, key=lambda x: x[1], reverse=True)
            best_match = sorted_list[0][0]
            count_loci.write(best_match + '\t' + locus + '\t' \
                             + str(len(loci_content[locus])) + '\n')

            if gene_overlap:
                for isoform in loci_content[locus]:
                    start, end, direction, blocksizes = isoform[0], isoform[1], isoform[2], isoform[3]
                    blockstarts, readstarts, name = isoform[4], isoform[5], isoform[6]
                    previous_blockend = start
                    previous_blockstart = start

                    for x in range(0, len(blocksizes), 1):
                        blockstart = int(blockstarts[x])
                        blocksize = int(blocksizes[x])
                        blockend = blockstart + blocksize
                        test = False

                        if blockstart - previous_blockend >= 50:
                             exonstart = previous_blockstart
                             exonend = previous_blockend
                             previous_blockstart = blockstart
                             previous_blockend = blockend
                             test = True

                        elif blockstart - previous_blockend < 50:
                             previous_blockend = blockend

                        if test and 'exon' in targets:
                             exon_dict = test_overlap(name, exonstart, exonend,
                                                      exon_dict, direction,
                                                      chromosome, new_exons)

                    exonstart = previous_blockstart
                    exonend = blockend
                    if 'exon' in targets:
                        exon_dict = test_overlap(name, exonstart, exonend,
                                                 exon_dict, direction,
                                                 chromosome, new_exons)

                    if not terminal_dict[chromosome].get(start):
                        if 'ends' in targets:
                            if direction == '+':
                                new_TSS.write(chromosome + '\t' + str(start) + '\t' \
                                              + str(start) + '\t' + name + '\t0\t' \
                                              + direction + '\n')
                            else:
                                new_polyA.write(chromosome + '\t' + str(start) \
                                                + '\t' + str(start) + '\t' + name \
                                                + '\t0\t' + direction + '\n')
                            for x in range(start-10, start+10, 1):
                                terminal_dict[chromosome][x] = 1

                    if not terminal_dict[chromosome].get(end):
                        if 'ends' in targets:
                            if direction == '-':
                                new_TSS.write(chromosome + '\t' + str(end) \
                                              + '\t' + str(end) + '\t' + name \
                                              + '\t0\t' + direction + '\n')
                            else:
                                new_polyA.write(chromosome + '\t' + str(end) \
                                                + '\t' + str(end) + '\t' + name \
                                                + '\t0\t' + direction + '\n')
                            for x in range(end-10, end+10, 1):
                                terminal_dict[chromosome][x] = 1

            else:
                if 'locus' in targets:
                    new_genes.write(chromosome + '\t' + str(locus_start) + '\t' \
                                    + str(locus_end) + '\t' + locus + ',')
                    for isoform in loci_content[locus]:
                        name = isoform[6]
                        new_genes.write(name + ',')
                    new_genes.write('\n')

def collect_all_chromosomes(gtf, isoforms):
    chromosome_set = set()
    for line in open(gtf):
        if line[0] != '#':
            a = line.strip().split('\t')
            chromosome_set.add(a[0])
    for line in open(isoforms):
        a = line.strip().split('\t')
        if len(a) > 16:
            chromosome_set.add(a[13])
    return sorted(list(chromosome_set))

def define_isoform_loci(isoform_file, target_chromosome):
    loci_dict, loci_content, counter = {}, {}, 1
    for line in open(isoform_file):
        a = line.strip().split('\t')
        chromosome = a[13]
        name = a[9]
        parts = name.split('_')
        splice = parts[1]
        if splice != '':
            if chromosome == target_chromosome:
                start = int(a[15])
                end = int(a[16])
                direction = a[8]
                name = a[9]

                begin, span = int(a[15]), int(a[16])
                blocksizes = a[18].split(',')[:-1]
                blockstarts = a[20].split(',')[:-1]
                readstarts = a[19].split(',')[:-1]
                overlap = False
                for base in range(begin, span, 1):
                    if not loci_dict.get(base):
                        pass
                    else:
                        locus = loci_dict[base]
                        loci_content[locus].append((start, end, direction,
                                                    blocksizes, blockstarts,
                                                    readstarts, name))
                        overlap = True
                        for base in range(begin, span, 1):
                            loci_dict[base] = locus
                        break
                if overlap == False:
                    locus = chromosome + '_' + str(counter)
                    for base in range(begin, span, 1):
                        loci_dict[base] = locus
                    loci_content[locus] = []
                    loci_content[locus].append((start, end, direction, blocksizes,
                                                blockstarts, readstarts, name))
                    counter += 1
    return loci_content

def count_exons(exon_file):
    counter = 0
    lengths = []
    for line in open(exon_file):
        counter += 1
        a = line.strip().split('\t')
        lengths.append(int(a[2]) - int(a[1]))
    return counter, lengths

def count_start_ends(tss_file):
    counter = 0
    for line in open(tss_file):
        counter += 1
    return counter

def count_new_loci(locus_file):
    counter, counter_spliced = 0, 0
    for line in open(locus_file):
        counter += 1
        a = line.strip().split('\t')
        isos = a[3].split(',')[1:-1]
        for iso in isos:
            if iso.split('_')[1] != '':
                counter_spliced += 1
    return counter, counter_spliced

def count_locus(locus_file):
    counter, isoforms = 0, []
    for line in open(locus_file):
        counter += 1
        isoforms.append(int(line.strip().split('\t')[3]))
    return counter, isoforms

def get_isoform_classification(classification):
    counter, dict1 = 0, {}
    type_set = set()
    for line in open(classification):
        counter += 1
        a = line.strip().split('\t')
        if counter > 1:
            dict1[a[0]] = a[5]
            type_set.add(a[5])

    for type1 in type_set:
        os.system('touch %s' %(path + '/' + type1 + '.psl'))
        os.system('touch %s' %(path + '/' + type1 + '.fasta'))

    return dict1

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

def split_reads(psl_file, class_dict, reads):
    out_psl_all = open(path + '/Isoforms_SQANTI_filtered.psl', 'w')
    out_fasta_all = open(path + '/Isoforms_SQANTI_filtered.fasta', 'w')

    for line in open(psl_file):
        a = line.strip().split('\t')
        name = a[9]
        if not class_dict.get(name):
            type1 = '-'
        else:
            type1 = class_dict[name]

        out_psl = open(path + '/' + type1 + '.psl', 'a')
        out_psl.write(line)
        out_psl.close()
        out_fasta = open(path + '/' + type1 + '.fasta', 'a')
        out_fasta.write('>' + name + '\n' + reads[name] + '\n')
        out_fasta.close()

        if type1 not in ['incomplete-splice_match']:
            out_psl_all.write(line)
            out_fasta_all.write('>' + name + '\n' + reads[name] + '\n')

def process_isoforms(chromosome_list, gtf_file, sqanti_folder):
    all_loci_verbose = open(sqanti_folder + 'R2C2_all_loci_verbose.txt', 'w')
    new_genes = open(sqanti_folder + 'R2C2_newly_identified_gene_loci.psl', 'w')
    new_exons = open(sqanti_folder + 'R2C2_newly_identified_exons.bed', 'w')
    new_TSS = open(sqanti_folder + 'R2C2_newly_identified_TSS.bed', 'w')
    new_polyA = open(sqanti_folder + 'R2C2_newly_identified_polyA.bed', 'w')
    all_loci = open(sqanti_folder + 'R2C2_all_loci.txt', 'w')
    count_loci = open(sqanti_folder + 'R2C2_count_loci.txt', 'w')

    combinations = [[[sqanti_folder + 'novel_not_in_catalog.psl',
                      sqanti_folder + '/novel_in_catalog.psl',
                      sqanti_folder + 'full-splice_match.psl'], ['ends']],
                    [[sqanti_folder + 'novel_not_in_catalog.psl'], ['exon']],
                    [[sqanti_folder + 'intergenic.psl'], ['locus']]]

    for combination in combinations:
            psl_files, targets = combination[0], combination[1]
            print(psl_files, targets)
    for chromosome in chromosome_list:
        print('gathering genome info')
        gene_dict, exon_dict, terminal_dict = read_genome_gtf(gtf_file, chromosome)
        print(chromosome)
        for combination in combinations:
            psl_files, targets = combination[0], combination[1]
            print(psl_files, targets)
            for psl_file in psl_files:
                print(psl_file)
                print('gathering read defined gene loci')
                loci_content = define_isoform_loci(psl_file, chromosome)
                print('evaluating isoforms')
                read_isoforms(loci_content, gene_dict, exon_dict, terminal_dict,
                              chromosome, targets, all_loci, new_genes, new_exons,
                              new_TSS, new_polyA, all_loci_verbose, count_loci)

    all_loci_verbose.close()
    new_genes.close()
    new_exons.close()
    new_TSS.close()
    new_polyA.close()
    all_loci.close()
    count_loci.close()

def print_summary_stats(sqanti_folder):
    counter, lengths = count_exons(sqanti_folder + 'R2C2_newly_identified_exons.bed')
    print(str(counter) + ' new exons of median length ' + str(np.average(lengths)))
    counter = count_start_ends(sqanti_folder + 'R2C2_newly_identified_TSS.bed')
    print(str(counter) + ' new TSSs')
    counter = count_start_ends(sqanti_folder + 'R2C2_newly_identified_polyA.bed')
    print(str(counter) + ' new polyAs')
    counter, counter_spliced = count_new_loci(sqanti_folder + 'R2C2_newly_identified_gene_loci.psl')
    counter_all, isoforms = count_locus(sqanti_folder+'R2C2_all_loci.txt')
    print(str(counter_all) + ' gene_loci with a total of ' + str(sum(isoforms)) \
          + ' isoforms.\n Of these loci ' + str(counter) \
          + ' do not overlap with known loci containing ' \
          + str(counter_spliced) + ' spliced isoforms')

def main():
    print('parsing SQANTI classification')
    class_dict = get_isoform_classification(classification)
    print('reading read fasta')
    reads = read_fasta(fasta_file)
    print('splitting reads based on SQANTI classification')
    split_reads(psl_file, class_dict, reads)
    print('collecting chromosomes')
    chromosome_list = collect_all_chromosomes(gtf_file, psl_file)
    print('finding new isoform features')
    process_isoforms(chromosome_list, gtf_file, path)
    print_summary_stats(path)

main()

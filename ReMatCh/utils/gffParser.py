#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import ntpath

version = '1.0'


def parse_id(filename):
    # get wanted feature IDs
    gff_ids = []
    with open(filename, 'r') as in_handle:
        for line in in_handle:
            line = line.strip()
            gff_ids.append(line)
    return gff_ids


def retrieve_seq_file(fasta_file, coord_file, extra_seq, filename, output_dir):
    # Parsing the sequence file, using the provided txt file containing the contig ID and positions to retrieve sequences.
    handle = open(fasta_file, "rU")
    records_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()

    seq_2_get = {}
    with open(coord_file, 'r') as sequeces2get:
        for line in sequeces2get:
            line = line.split(',')
            coords = (int(line[-2]), int(line[-1]))
            contig_id = line[0]
            if contig_id in list(seq_2_get.keys()):
                seq_2_get[contig_id].append(coords)
            else:
                seq_2_get[contig_id] = [coords]

    with open(output_dir + '/' + filename + '.fasta', 'w') as output_handle:
        fails = 0
        successes = 0
        records = []
        for contig, listCoords in list(seq_2_get.items()):
            contig_seq = records_dict[contig].seq
            for coord in listCoords:
                coord1 = coord[0] - extra_seq
                coord2 = coord[1] + extra_seq
                if coord1 < 0 or coord2 > len(contig_seq):
                    fail_log = open(output_dir + '/' + filename + '_fails.txt', 'a')
                    fail_log.write(contig + ',' + str(coord[0]) + ',' + str(coord[1]) + '\n')
                    fail_log.close()
                    fails += 1
                else:
                    geneseq = str(contig_seq[coord1:coord2])
                    record = SeqRecord(Seq(geneseq), id=str(str(contig) + '#' + str(coord1) + '_' + str(coord2)),
                                       description='')
                    records.append(record)
                    successes += 1
        SeqIO.write(records, output_handle, "fasta")

    print('Retrived %s features successfully from %s with %s bp as extra'
          ' sequence.' % (str(successes), filename, str(extra_seq)))
    if fails > 0:
        print('%s featrued failed to retrieve. Check %s_fails.txt file.' % (str(fails), filename))


def retrieve_seq(fasta_file, gff_features, extra_seq, filename, output_dir):
    # parsing the sequence file into a SeqIO dictionary. one contig per entry
    handle = open(fasta_file, "rU")
    records_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()

    with open(output_dir + '/' + filename + '.fasta', 'w') as output_handle:
        fails = 0
        successes = 0
        records = []
        for locus, location in list(gff_features.items()):
            # print locus
            contig_seq = records_dict[location[0]].seq
            coord1 = location[1] - extra_seq
            coord2 = location[2] + extra_seq
            if coord1 < 0 or coord2 > len(contig_seq):
                fail_log = open(output_dir + '/' + filename + '_fails.txt', 'a')
                fail_log.write(locus + '\n')
                fail_log.close()
                fails += 1
            else:
                geneseq = str(contig_seq[coord1:coord2])
                if location[3] == '-':
                    seq = Seq(geneseq)
                    geneseq = str(seq.reverse_complement())
                record = SeqRecord(Seq(geneseq),
                                   id=str(locus + '-' + str(location[0]) + '#' + str(location[1]) + '_' +
                                          str(location[2])),
                                   description='')
                records.append(record)
                successes += 1
        SeqIO.write(records, output_handle, "fasta")
    print('Retrived %s features successfully from %s with %s bp as extra'
          ' sequence.' % (str(successes), filename, str(extra_seq)))
    if fails > 0:
        print('%s featrued failed to retrieve. Check %s_fails.txt file.' % (str(fails), filename))


def parse_features(temp_gff):
    # parsing the feature file into a dictionary
    gff_features = {}

    with open(temp_gff, 'r') as temp_genes:
        for line in temp_genes:
            line = line.split('\t')
            if "CDS" in line[2]:
                id = line[-1].split(';')
                locus_id = str(id[0].split('=')[1])
                contig = line[0]
                begining = int(line[3]) - 1  # to get the full sequence
                end = int(line[4])
                strand = line[6]
                location = [contig, begining, end, strand]
                gff_features[locus_id] = location
    return gff_features


def gff_parser(gff_file, extra_seq=0, output_dir='.', keep_temporary_files=False, ids=None, coord_file=None):
    filename = ntpath.basename(gff_file).replace('.gff', '')

    # cleaning temp files if they exist
    if os.path.isfile(output_dir + '/' + filename + '_features.gff'):
        os.remove(output_dir + '/' + filename + '_features.gff')
    if os.path.isfile(output_dir + '/' + filename + '_sequence.fasta'):
        os.remove(output_dir + '/' + filename + '_sequence.fasta')

    # cleaning fails file if it exists
    if os.path.isfile(output_dir + '/' + filename + '_fails.txt'):
        os.remove(output_dir + '/' + filename + '_fails.txt')

    if coord_file is None:

        if ids is not None:
            select_ids = parse_id(ids)
        else:
            select_ids = None

        # separating the gff into 2 different files: one with the features and another with the conting sequences
        with open(gff_file, 'r') as in_handle, open(output_dir + '/' + filename + '_features.gff', 'a') as temp_genes, \
                open(output_dir + '/' + filename + '_sequence.fasta', 'a') as temp_contigs:
            for line in in_handle:
                if not line.startswith('##'):
                    if '\t' in line:
                        if select_ids is not None:
                            items = line.split('\t')
                            id = items[-1].split(';')[0]
                            id = id.split('=')[1]
                            if id in select_ids:
                                temp_genes.write(line)
                        else:
                            temp_genes.write(line)
                    else:
                        temp_contigs.write(line)

        gff_files = parse_features(output_dir + '/' + filename + '_features.gff')

        retrieve_seq(output_dir + '/' + filename + '_sequence.fasta', gff_files, extra_seq, filename, output_dir)

    else:
        with open(gff_file, 'r') as in_handle, \
                open(output_dir + '/' + filename + '_sequence.fasta', 'a') as temp_contigs:
            for line in in_handle:
                if not line.startswith('##'):
                    if '\t' in line:
                        pass
                    else:
                        temp_contigs.write(line)

        retrieve_seq_file(output_dir + '/' + filename + '_sequence.fasta', coord_file, extra_seq, filename, output_dir)

    # removing temp files
    if not keep_temporary_files:
        try:
            os.remove(output_dir + '/' + filename + '_features.gff')
        except:
            pass
        os.remove(output_dir + '/' + filename + '_sequence.fasta')


def main():
    parser = argparse.ArgumentParser(prog='gffParser.py', description='GFF3 parser for feature sequence retrival.',
                                     epilog='by C I Mendes (cimendes@medicina.ulisboa.pt)')
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser.add_argument('-i', '--input',
                        help='GFF3 file to parse, containing both sequences and annotations (like the one obtained from'
                             ' PROKKA).', type=argparse.FileType('r'), required=True)
    parser.add_argument('-x', '--extraSeq', help='Extra sequence to retrieve per feature in gff.', default=0, type=int,
                        required=False)
    parser.add_argument('-k', '--keepTemporaryFiles', help='Keep temporary gff(without sequence) and fasta files.',
                        action='store_true')
    parser.add_argument('-o', '--outputDir', help='Path to where the output is to be saved.', default='.',
                        required=False)

    parser_optional_selected_regions_exclusive = parser.add_mutually_exclusive_group()
    parser_optional_selected_regions_exclusive.add_argument('-s', '--select',
                                                            help='txt file with the IDs of interest, one per line',
                                                            type=argparse.FileType('r'), required=False)
    parser_optional_selected_regions_exclusive.add_argument('-f', '--fromFile',
                                                            help='Sequence coordinates to be retrieved. Requires contig'
                                                                 ' ID and coords (contig,strart,end) in a csv file. One'
                                                                 ' per line.', type=argparse.FileType('r'),
                                                            required=False)

    args = parser.parse_args()

    args.outputDir = os.path.abspath(args.outputDir)
    if not os.path.isdir(args.outputDir):
        os.makedirs(args.outputDir)

    gff_parser(os.path.abspath(args.input.name), args.extraSeq, os.path.abspath(args.outputDir),
               args.keepTemporaryFiles,
               os.path.abspath(args.select.name) if args.select is not None else None,
               os.path.abspath(args.fromFile.name) if args.fromFile is not None else None)


if __name__ == "__main__":
    main()

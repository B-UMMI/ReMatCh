#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
strip_alignment.py - Strip alignment positions containing gaps,
missing data and invariable positions
<https://github.com/B-UMMI/ReMatCh/>

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: October 15, 2018

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from Bio import SeqIO
import os
import argparse
import sys


version = '0.2'


def get_sequences(infile):
    print('Getting sequences')
    sequences_seq_io = list(SeqIO.parse(infile, 'fasta'))

    sequence_length = None
    sequences_dict = {}
    all_executed_printed = False
    for x, sequence in enumerate(sequences_seq_io):
        if sequence_length is None:
            sequence_length = len(sequence.seq)
        if sequence_length != len(sequence.seq):
            sys.exit('Sequences with different length!')
        sequences_dict[sequence.id] = list(sequence.seq)

        if (x + 1) % 10 == 0:
            print('\n' + str(round((float(x + 1) / len(sequences_seq_io)) * 100, 2)) + '% of sequences already processed (getting sequences)')
            if x + 1 == len(sequences_seq_io):
                all_executed_printed = True
    if not all_executed_printed:
        print('\n' + str(round((float(x + 1) / len(sequences_seq_io)) * 100, 2)) + '% of sequences already processed (getting sequences)')

    return sequences_dict, sequence_length


def positions_type(sequences_dict, sequence_length, not_gaps, not_missing, not_invariable):
    print('Determining positions type')
    positions_2_keep = []
    invariable = []
    missing = []
    gaps = []
    gaps_missing = 0
    all_executed_printed = False
    for i in range(0, sequence_length):
        data = []
        for sample in sequences_dict:
            data.append(sequences_dict[sample][i])
        possibilities = set(data)
        if len(possibilities) == 1:
            invariable.append(i)
        if len(possibilities.intersection(set(['N']))) > 0:
            missing.append(i)
        if len(possibilities.intersection(set(['-']))) > 0:
            gaps.append(i)
        if len(possibilities.intersection(set(['N', '-']))) > 0:
            gaps_missing += 1
        if len(possibilities) > 1 and len(possibilities.intersection(set(['N', '-']))) == 0:
            positions_2_keep.append(i)

        if (i + 1) % 10000 == 0:
            print('\n' + str(round((float(i + 1) / sequence_length) * 100, 2)) + '% of positions already'
                                                                                 ' processed (determining positions'
                                                                                 ' type)')
            if i + 1 == len(sequences_dict):
                all_executed_printed = True
        if not all_executed_printed:
            print('\n' + str(round((float(i + 1) / sequence_length) * 100, 2)) + '% of positions already'
                                                                                 ' processed (determining positions'
                                                                                 ' type)')

    print('Positions to keep (no matter): ' + str(len(positions_2_keep)))
    print('Invariable sites: ' + str(len(invariable)))
    print('Positions with missing data ("N"): ' + str(len(missing)))
    print('Positions with GAPs ("-"): ' + str(len(gaps)))
    print('Positions with GAPs or missing data: ' + str(gaps_missing))

    if not_gaps:
        positions_2_keep.extend(gaps)
    if not_missing:
        positions_2_keep.extend(missing)
    if not_invariable:
        positions_2_keep.extend(invariable)

    positions_2_keep = sorted(set(positions_2_keep))

    print('Positions to keep (final): ' + str(len(positions_2_keep)))

    return positions_2_keep


def chunkstring(string, length):
    return (string[0 + i:length + i] for i in range(0, len(string), length))


def write_fasta(sequences_dict, positions_2_keep, outfile):
    print('Writing stripped sequences')
    all_executed_printed = False
    with open(outfile, 'wt') as writer:
        for x, sample in enumerate(sequences_dict):
            writer.write('>' + sample + '\n')
            fasta_sequence_lines = chunkstring(''.join([sequences_dict[sample][i] for i in positions_2_keep]), 80)
            for line in fasta_sequence_lines:
                writer.write(line + '\n')

            if (x + 1) % 100 == 0:
                print('\n' + str(round((float(x + 1) / len(sequences_dict)) * 100, 2)) + '% of sequences already'
                                                                                         ' processed (writing stripped'
                                                                                         ' sequences)')
                if x + 1 == len(sequences_dict):
                    all_executed_printed = True
        if not all_executed_printed:
            print('\n' + str(round((float(x + 1) / len(sequences_dict)) * 100, 2)) + '% of sequences already'
                                                                                     ' processed (writing stripped'
                                                                                     ' sequences)')


def strip_alignment(args):
    outdir = os.path.dirname(os.path.abspath(args.outfile))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    outfile = os.path.abspath(args.outfile)

    infile = os.path.abspath(args.infile.name)

    sequences_dict, sequence_length = get_sequences(infile)
    positions_2_keep = positions_type(sequences_dict, sequence_length, args.notGAPs, args.notMissing,
                                      args.notInvariable)
    write_fasta(sequences_dict, positions_2_keep, outfile)


def main():
    parser = argparse.ArgumentParser(prog='strip_alignment.py',
                                     description='Strip alignment positions containing gaps, missing data and'
                                                 ' invariable positions',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-i', '--infile', type=argparse.FileType('r'),
                                 metavar='/path/to/aligned/input/file.fasta', help='Path to the aligned fasta file',
                                 required=True)
    parser_required.add_argument('-o', '--outfile', type=str, metavar='/path/to/stripped/output/file.fasta',
                                 help='Stripped output fasta file', required=True, default='alignment_stripped.fasta')

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('--notGAPs', action='store_true', help='Not strip positions with GAPs')
    parser_optional_general.add_argument('--notMissing', action='store_true',
                                         help='Not strip positions with missing data')
    parser_optional_general.add_argument('--notInvariable', action='store_true', help='Not strip invariable sites')

    args = parser.parse_args()

    strip_alignment(args)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
combine_alignment_consensus.py - Combine the alignment consensus
sequences from ReMatCh first run by reference sequences into single
files
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

import os
import argparse
import time
import sys

version = '0.2'


def concatenate_files(input_files_list, outdir):
    all_executed_printed = False
    for x, input_file in enumerate(input_files_list):
        sample = os.path.basename(input_file).rsplit('.', 2)[0]
        with open(input_file, 'rtU') as reader:
            writer = None
            for line in reader:
                line = line.rstrip('\r\n')
                if line.startswith('>'):
                    file_output = os.path.join(outdir, line[1:] + '.fasta')
                    if writer is not None:
                        writer.flush()
                        writer.close()
                    if os.path.isfile(file_output):
                        writer = open(file_output, 'at')
                    else:
                        writer = open(file_output, 'wt')
                    writer.write('>' + sample + '\n')
                else:
                    if len(line) > 0:
                        writer.write(line + '\n')
            writer.flush()
            writer.close()

        if (x + 1) % 100 == 0:
            print('\n' + str(round((float(x + 1) / len(input_files_list)) * 100, 2)) + '% of IDs already processed')
            all_executed_printed = True
    if not all_executed_printed:
        print('\n' + str(round((float(x + 1) / len(input_files_list)) * 100, 2)) + '% of IDs already processed')


def combine_alignment_consensus(args):
    outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    outdir = os.path.join(outdir, 'combine_alignment_consensus_' + time.strftime("%Y%m%d-%H%M%S"), '')
    os.makedirs(outdir)

    workdir = os.path.abspath(args.workdir)

    alignment_files = []
    directories = [d for d in os.listdir(workdir) if
                   not d.startswith('.') and
                   os.path.isdir(os.path.join(workdir, d, ''))]
    for sample_dir in directories:
        sample_dir_path = os.path.join(workdir, sample_dir, '')
        files = [f for f in os.listdir(sample_dir_path) if
                 not f.startswith('.') and
                 os.path.isfile(os.path.join(sample_dir_path, f))]
        for file_found in files:
            if file_found.endswith('.alignment.fasta'):
                file_found_path = os.path.join(sample_dir_path, file_found)
                alignment_files.append(file_found_path)

    if len(alignment_files) > 0:
        concatenate_files(alignment_files, outdir)
    else:
        sys.exit('No ReMatCh alignment.fasta files were found!')


def main():
    parser = argparse.ArgumentParser(prog='combine_alignment_consensus.py',
                                     description='Combine the alignment consensus sequences from ReMatCh first run by'
                                                 ' reference sequences into single'
                                                 ' files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-w', '--workdir', type=str, metavar='/path/to/rematch/working/directory/',
                                 help='Path to the directory where ReMatCh was running', required=True)

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                         help='Path to the directory where the combined sequence files will stored',
                                         required=False, default='.')

    args = parser.parse_args()

    combine_alignment_consensus(args)


if __name__ == "__main__":
    main()

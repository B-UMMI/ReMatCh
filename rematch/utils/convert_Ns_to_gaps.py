#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
convert_Ns_to_gaps.py - Convert the Ns into gaps
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


version = '0.2'


def conversion(infile, outfile):
    last_printed = 0
    counter = 1
    with open(infile, 'rtU') as reader:
        with open(outfile, 'wt') as writer:
            for line in reader:
                line = line.rstrip('\r\n')
                if line.startswith('>'):
                    writer.write(line + '\n')
                    if counter % 10 == 0:
                        print('\n' + str(counter) + ' sequences already processed')
                        last_printed = counter
                    counter += 1
                else:
                    if len(line) > 0:
                        line = line.replace('N', '-')
                        writer.write(line + '\n')
    if last_printed < counter:
        print('\n' + str(counter - 1) + ' sequences already processed')


def convert_n_2_gaps(args):
    outdir = os.path.dirname(os.path.abspath(args.outfile))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    outfile = os.path.abspath(args.outfile)

    infile = os.path.abspath(args.infile.name)

    conversion(infile, outfile)


def main():
    parser = argparse.ArgumentParser(prog='convert_Ns_to_gaps.py', description='Convert the Ns into gaps',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-i', '--infile', type=argparse.FileType('r'), metavar='/path/to/input/file.fasta',
                                 help='Path to the fasta file', required=True)
    parser_required.add_argument('-o', '--outfile', type=str, metavar='/path/to/converted/output/file.fasta',
                                 help='Converted output fasta file', required=True,
                                 default='converted_Ns_to_gaps.fasta')

    args = parser.parse_args()

    convert_n_2_gaps(args)


if __name__ == "__main__":
    main()

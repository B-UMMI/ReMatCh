#!/usr/bin/env python3

# -*- coding: utf-8 -*-

'''
Adapted from:
https://github.com/mickaelsilva/pythonscripts/blob/master/SeqOfWeb/SeqFromWebTaxon.py
mickaelsilva
'''

import sys
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import time
import argparse
import os


def run_seq_from_web_taxon(taxonname, outputfile, getmachine, getOmicsDataType, getLibraryType, print_True):
    print('\n' + 'Searching RunIDs for ' + taxonname)

    taxonname = urllib.parse.quote(taxonname)
    url = "http://www.ebi.ac.uk/ena/data/view/Taxon%3A" + taxonname + "&display=xml"
    try:
        content = urllib.request.urlopen(url)
        xml = content.read()
        tree = ET.fromstring(xml)
        taxonid = ''
    except:
        print("Ooops!There might be a problem with the ena service, try later or check if the xml is well formated"
              " at " + url)
        raise
    for child in tree:
        taxonid = child.get('taxId')
    if (taxonid):
        print("\n" + "Taxon ID found: " + taxonid)
        url = "http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_tree%28" + \
              taxonid + \
              "%29%22&result=read_run&display=xml"

        content = urllib.request.urlopen(url)
        xml = content.read()
        tree = ET.fromstring(xml)

        runid = ''
        n = 0
        with open(outputfile, "wt") as f:
            f.write('#' + str(time.strftime("%d/%m/%Y")) + "\n")
            model = ''
            prjid = ''
            length_line = 0
            omics = ''
            libraryType = ''
            for child in tree:
                runid = child.get('accession')

                n += 1

                if getmachine is True or getOmicsDataType is True or getLibraryType is True:
                    for child2 in child:
                        if child2.tag == 'EXPERIMENT_REF':
                            expid = child2.get('accession')
                            url2 = "http://www.ebi.ac.uk/ena/data/view/" + expid + "&display=xml"
                            content = urllib.request.urlopen(url2)
                            xml = content.read()
                            tree2 = ET.fromstring(xml)
                            try:
                                for child3 in tree2:
                                    for child4 in child3:
                                        if child4.tag == 'PLATFORM':
                                            for child5 in child4:
                                                for child6 in child5:
                                                    if child6.tag == 'INSTRUMENT_MODEL':
                                                        model = child6.text
                                        elif child4.tag == 'STUDY_REF':
                                            prjid = child4.get('accession')
                                        elif child4.tag == 'DESIGN':
                                            if getOmicsDataType is True or getLibraryType is True:
                                                for child5 in child4:
                                                    if child5.tag == 'LIBRARY_DESCRIPTOR':
                                                        for child6 in child5:
                                                            if child6.tag == 'LIBRARY_SOURCE' and getOmicsDataType is True:
                                                                omics = child6.text
                                                            elif child6.tag == 'LIBRARY_LAYOUT' and getLibraryType is True:
                                                                libraryType = child6[0].tag
                            except:
                                model = 'not found'
                                omics = 'not found'
                                libraryType = 'not found'
                    f.write(str(runid) + "\t" + model + "\t" + prjid + "\t" + omics + "\t" + libraryType + "\n")
                    if print_True:
                        line = "run acession %s sequenced on %s from project %s for %s %s end" \
                               " data" % (runid, model, prjid, omics, libraryType)
                        if length_line < len(line):
                            length_line = len(line)
                        sys.stderr.write("\r" + line + str(' ' * (length_line - len(line))))
                        sys.stderr.flush()
                else:
                    f.write(str(runid) + '\t' * 4 + "\n")
                    if print_True:
                        line = "run acession %s" % (runid, prjid)
                        if length_line < len(line):
                            length_line = len(line)
                        sys.stderr.write("\r" + line + str(' ' * (length_line - len(line))))
                        sys.stderr.flush()
        print("\n")
        print("\n"
              "found %s run id's" % n)

    else:
        print("taxon name does not exist")


def main():
    parser = argparse.ArgumentParser(description="This program gets a list of sequencing runs and machine were the"
                                                 " sequencing was performed, given a taxon name accepted by the"
                                                 " European nucleotide Archive")
    parser.add_argument('-i', nargs=1, type=str, help='taxon name', metavar='"Streptococcus agalactiae"', required=True)
    parser.add_argument('-o', nargs=1, type=str, help='output file name', required=True)
    parser.add_argument('-g', help='True to include sequencing machine in the output', action='store_true',
                        required=False)
    parser.add_argument('--getOmicsDataType', help='Informs the programme to include OMICS data type'
                                                   ' (examples: GENOMIC / TRANSCRIPTOMIC / SYNTHETIC) in the output',
                        action='store_true')
    parser.add_argument('--getLibraryType', help='Informs the programme to include library type'
                                                 ' (examples: PAIRED / SINGLE) in the output', action='store_true')

    args = parser.parse_args()

    getmachine = args.g
    taxonname = args.i[0]

    outdir = os.path.dirname(os.path.abspath(args.o[0]))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outputfile = os.path.abspath(args.o[0])

    getOmicsDataType = args.getOmicsDataType
    getLibraryType = args.getLibraryType

    run_seq_from_web_taxon(taxonname, outputfile, getmachine, getOmicsDataType, getLibraryType, True)


if __name__ == "__main__":
    main()

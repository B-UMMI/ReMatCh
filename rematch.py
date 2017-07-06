#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
rematch.py - Reads mapping against target sequences, checking mapping
and consensus sequences production
<https://github.com/B-UMMI/ReMatCh/>

Copyright (C) 2017 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: April 12, 2017

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
import sys
import time
import argparse
import modules.utils as utils
import modules.seqFromWebTaxon as seqFromWebTaxon
import modules.download as download
import modules.rematch_module as rematch_module
import modules.checkMLST as checkMLST


version = '3.2'


def searchFastqFiles(directory):
    filesExtensions = ['.fastq.gz', '.fq.gz']
    pairEnd_filesSeparation = [['_R1_001.f', '_R2_001.f'], ['_1.f', '_2.f']]

    listIDs = {}
    directories = [d for d in os.listdir(directory) if not d.startswith('.') and os.path.isdir(os.path.join(directory, d, ''))]
    for directory_found in directories:
        if directory_found != 'pubmlst':
            directory_path = os.path.join(directory, directory_found, '')

            fastqFound = []
            files = [f for f in os.listdir(directory_path) if not f.startswith('.') and os.path.isfile(os.path.join(directory_path, f))]
            for file_found in files:
                if file_found.endswith(tuple(filesExtensions)):
                    fastqFound.append(file_found)

            if len(fastqFound) == 1:
                listIDs[directory_found] = [os.path.join(directory_path, f) for f in fastqFound]
            elif len(fastqFound) >= 2:
                file_pair = []

                # Search pairs
                for PE_separation in pairEnd_filesSeparation:
                    for fastq in fastqFound:
                        if PE_separation[0] in fastq or PE_separation[1] in fastq:
                            file_pair.append(fastq)

                    if len(file_pair) == 2:
                        break
                    else:
                        file_pair = []

                # Search single
                if len(file_pair) == 0:
                    for PE_separation in pairEnd_filesSeparation:
                        for fastq in fastqFound:
                            if PE_separation[0] not in fastq or PE_separation[1] not in fastq:
                                file_pair.append(fastq)

                    if len(file_pair) >= 1:
                        file_pair = file_pair[0]

                if len(file_pair) >= 1:
                    listIDs[directory_found] = [os.path.join(directory_path, f) for f in file_pair]

    return listIDs


def getListIDs_fromFile(fileListIDs):
    list_ids = []

    with open(fileListIDs, 'rtU') as lines:
        for line in lines:
            line = line.splitlines()[0]
            if len(line) > 0:
                list_ids.append(line)

    if len(list_ids) == 0:
        sys.exit('No runIDs were found in ' + fileListIDs)

    return list_ids


def getTaxonRunIDs(taxon_name, outputfile):
    seqFromWebTaxon.runSeqFromWebTaxon(taxon_name, outputfile, True, True, True, False)

    runIDs = []
    with open(outputfile, 'rtU') as reader:
        for line in reader:
            line = line.splitlines()[0]
            if len(line) > 0:
                if not line.startswith('#'):
                    line = line.split('\t')
                    runIDs.append(line[0])

    return runIDs


def getListIDs(workdir, fileListIDs, taxon_name):
    searched_fastq_files = False
    listIDs = []
    if fileListIDs is None and taxon_name is None:
        listIDs = searchFastqFiles(workdir)
        searched_fastq_files = True
    elif fileListIDs is not None:
        listIDs = getListIDs_fromFile(os.path.abspath(fileListIDs))
    elif taxon_name is not None and fileListIDs is None:
        listIDs = getTaxonRunIDs(taxon_name, os.path.join(workdir, 'IDs_list.seqFromWebTaxon.tab'))

    if len(listIDs) == 0:
        sys.exit('No IDs were found')
    return listIDs, searched_fastq_files


def format_gene_info(gene_specific_info, minimum_gene_coverage, minimum_gene_identity, reported_data_type):
    info = None
    if gene_specific_info['gene_coverage'] >= minimum_gene_coverage and gene_specific_info['gene_identity'] >= minimum_gene_identity:
        if gene_specific_info['gene_number_positions_multiple_alleles'] == 0:
            info = str(gene_specific_info[reported_data_type])
        else:
            info = 'multiAlleles_' + str(gene_specific_info[reported_data_type])
    else:
        info = 'absent_' + str(gene_specific_info[reported_data_type])

    return info


def write_data_by_gene(gene_list_reference, minimum_gene_coverage, sample, data_by_gene, outdir, time_str, run_times, minimum_gene_identity, reported_data_type):
    combined_report = os.path.join(outdir, 'combined_report.data_by_gene.' + run_times + '.' + reported_data_type + '.' + time_str + '.tab')

    if reported_data_type == 'coverage_depth':
        reported_data_type = 'gene_mean_read_coverage'
    elif reported_data_type == 'sequence_coverage':
        reported_data_type = 'gene_coverage'

    combined_report_exist = os.path.isfile(combined_report)
    with open(combined_report, 'at') as writer:
        seq_list = gene_list_reference.keys()
        if not combined_report_exist:
            writer.write('#sample' + '\t' + '\t'.join([gene_list_reference[seq] for seq in seq_list]) + '\n')

        results = {}
        headers = []
        for i in data_by_gene:
            results[data_by_gene[i]['header']] = format_gene_info(data_by_gene[i], minimum_gene_coverage, minimum_gene_identity, reported_data_type)
            headers.append(data_by_gene[i]['header'])

        if len(headers) != gene_list_reference:
            for gene in gene_list_reference:
                if gene not in headers:
                    results[gene] = 'NA'

        writer.write(sample + '\t' + '\t'.join([results[seq] for seq in seq_list]) + '\n')


def write_sample_report(sample, outdir, time_str, fileSize, run_successfully_fastq, run_successfully_rematch_first, run_successfully_rematch_second, time_taken_fastq, time_taken_rematch_first, time_taken_rematch_second, time_taken_sample, sequencingInformation, sample_data_general_first, sample_data_general_second, fastq_used):
    sample_report = os.path.join(outdir, 'sample_report.' + time_str + '.tab')
    report_exist = os.path.isfile(sample_report)

    header_general = ['sample', 'sample_run_successfully', 'sample_run_time', 'files_size', 'download_run_successfully', 'download_run_time', 'rematch_run_successfully_first', 'rematch_run_time_first', 'rematch_run_successfully_second', 'rematch_run_time_second']
    header_data_general = ['number_absent_genes', 'number_genes_multiple_alleles', 'mean_sample_coverage']
    header_sequencing = ['run_accession', 'instrument_platform', 'instrument_model', 'library_layout', 'library_source', 'extra_run_accession', 'nominal_length', 'read_count', 'base_count', 'date_download']

    with open(sample_report, 'at') as writer:
        if not report_exist:
            writer.write('#' + '\t'.join(header_general) + '\t' + '_first\t'.join(header_data_general) + '_first\t' + '_second\t'.join(header_data_general) + '_second\t' + '\t'.join(header_sequencing) + '\t' + 'fastq_used' + '\n')

        writer.write('\t'.join([sample, str(all([run_successfully_fastq is not False, run_successfully_rematch_first is not False, run_successfully_rematch_second is not False])), str(time_taken_sample), str(fileSize), str(run_successfully_fastq), str(time_taken_fastq), str(run_successfully_rematch_first), str(time_taken_rematch_first), str(run_successfully_rematch_second), str(time_taken_rematch_second)]) + '\t' + '\t'.join([str(sample_data_general_first[i]) for i in header_data_general]) + '\t' + '\t'.join([str(sample_data_general_second[i]) for i in header_data_general]) + '\t' + '\t'.join([str(sequencingInformation[i]) for i in header_sequencing]) + '\t' + ','.join(fastq_used) + '\n')


def concatenate_extraSeq_2_consensus(consensus_sequence, reference_sequence, extraSeq_length, outdir):
    reference_dict, ignore, ignore = rematch_module.get_sequence_information(reference_sequence, extraSeq_length)
    consensus_dict, genes, ignore = rematch_module.get_sequence_information(consensus_sequence, 0)
    for k, values_consensus in consensus_dict.items():
        for values_reference in reference_dict.values():
            if values_reference['header'] == values_consensus['header']:
                if extraSeq_length <= len(values_reference['sequence']):
                    right_extra_seq = '' if extraSeq_length == 0 else values_reference['sequence'][-extraSeq_length:]
                    consensus_dict[k]['sequence'] = values_reference['sequence'][:extraSeq_length] + consensus_dict[k]['sequence'] + right_extra_seq
                    consensus_dict[k]['length'] += extraSeq_length + len(right_extra_seq)

    consensus_concatenated = os.path.join(outdir, 'consensus_concatenated_extraSeq.fasta')
    with open(consensus_concatenated, 'wt') as writer:
        for i in consensus_dict:
            writer.write('>' + consensus_dict[i]['header'] + '\n')
            fasta_sequence_lines = rematch_module.chunkstring(consensus_dict[i]['sequence'], 80)
            for line in fasta_sequence_lines:
                writer.write(line + '\n')

    return consensus_concatenated, genes, consensus_dict


def clean_headers_reference_file(reference_file, outdir, extraSeq):
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", ":"]
    print 'Checking if reference sequences contain ' + str(problematic_characters) + '\n'
    headers_changed = False
    new_reference_file = reference_file
    sequences, genes, headers_changed = rematch_module.get_sequence_information(reference_file, extraSeq)
    if headers_changed:
        print 'At least one of the those characters was found. Replacing those with _' + '\n'
        new_reference_file = os.path.join(outdir, os.path.splitext(os.path.basename(reference_file))[0] + '.headers_renamed.fasta')
        with open(new_reference_file, 'wt') as writer:
            for i in sequences:
                writer.write('>' + sequences[i]['header'] + '\n')
                fasta_sequence_lines = rematch_module.chunkstring(sequences[i]['sequence'], 80)
                for line in fasta_sequence_lines:
                    writer.write(line + '\n')
    return new_reference_file, genes, sequences


def write_mlst_report(sample, run_times, consensus_type, st, alleles_profile, lociOrder, outdir, time_str):
    mlst_report = os.path.join(outdir, 'mlst_report.' + time_str + '.tab')
    mlst_report_exist = os.path.isfile(mlst_report)
    with open(mlst_report, 'at') as writer:
        if not mlst_report_exist:
            writer.write('\t'.join(['#sample', 'ReMatCh_run', 'consensus_type', 'ST'] + lociOrder) + '\n')
        writer.write('\t'.join([sample, run_times, consensus_type, str(st)] + alleles_profile.split(',')) + '\n')


def run_get_st(sample, mlst_dicts, consensus_sequences, mlstConsensus, run_times, outdir, time_str):
    if mlstConsensus == 'all':
        for consensus_type in consensus_sequences:
            print 'Searching MLST for ' + consensus_type + ' consensus'
            st, alleles_profile = checkMLST.getST(mlst_dicts, consensus_sequences[consensus_type])
            write_mlst_report(sample, run_times, consensus_type, st, alleles_profile, mlst_dicts[2], outdir, time_str)
            print 'ST found: ' + str(st) + ' (' + alleles_profile + ')'
    else:
        st, alleles_profile = checkMLST.getST(mlst_dicts, consensus_sequences[mlstConsensus])
        write_mlst_report(sample, run_times, mlstConsensus, st, alleles_profile, mlst_dicts[2], outdir, time_str)
        print 'ST found for ' + mlstConsensus + ' consensus: ' + str(st) + ' (' + alleles_profile + ')'


def runRematch(args):
    workdir = os.path.abspath(args.workdir)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    asperaKey = os.path.abspath(args.asperaKey.name) if args.asperaKey is not None else None

    # Start logger
    logfile, time_str = utils.start_logger(workdir)

    # Get general information
    script_path = utils.general_information(logfile, version, workdir, time_str, args.doNotUseProvidedSoftware, asperaKey, args.downloadCramBam)

    # Set listIDs
    listIDs, searched_fastq_files = getListIDs(workdir, args.listIDs.name if args.listIDs is not None else None, args.taxon)

    if args.mlst is not None:
        time_taken_PubMLST, mlst_dicts, mlst_sequences = checkMLST.downloadPubMLSTxml(args.mlst, args.mlstSchemaNumber, workdir)
        args.softClip_recodeRun = 'first'
        args.conservedSeq = False

    if args.reference is None:
        reference_file = checkMLST.check_existing_schema(args.mlst, args.mlstSchemaNumber, script_path)
        args.extraSeq = 200
        if reference_file is None:
            print 'It was not found provided MLST scheme sequences for ' + args.mlst
            print 'Trying to obtain reference MLST sequences from PubMLST'
            if len(mlst_sequences) > 0:
                reference_file = checkMLST.write_mlst_reference(args.mlst, mlst_sequences, workdir, time_str)
                args.extraSeq = 0
            else:
                sys.exit('It was not possible to download MLST sequences from PubMLST!')
        else:
            print 'Using provided scheme as referece: ' + reference_file
    else:
        reference_file = os.path.abspath(args.reference.name)

    # Run ReMatCh for each sample
    print '\n' + 'STARTING ReMatCh' + '\n'

    # Clean sequences headers
    reference_file, gene_list_reference, reference_dict = clean_headers_reference_file(reference_file, workdir, args.extraSeq)

    if args.mlst is not None:
        problem_genes = False
        for header in mlst_sequences:
            if header not in gene_list_reference:
                print 'MLST gene {header} not found between reference sequences'.format(header=header)
                problem_genes = True
        if problem_genes:
            sys.exit('Missing MLST genes from reference sequences (at least sequences names do not match)!')

    if len(gene_list_reference) == 0:
        sys.exit('No sequences left')

    # To use in combined report

    number_samples_successfully = 0
    for sample in listIDs:
        sample_start_time = time.time()
        print '\n\n' + 'Sample ID: ' + sample

        # Create sample outdir
        sample_outdir = os.path.join(workdir, sample, '')
        if not os.path.isdir(sample_outdir):
            os.mkdir(sample_outdir)

        run_successfully_fastq = None
        time_taken_fastq = 0
        sequencingInformation = {'run_accession': None, 'instrument_platform': None, 'instrument_model': None, 'library_layout': None, 'library_source': None, 'extra_run_accession': None, 'nominal_length': None, 'read_count': None, 'base_count': None, 'date_download': None}
        if not searched_fastq_files:
            # Download Files
            time_taken_fastq, run_successfully_fastq, fastq_files, sequencingInformation = download.runDownload(sample, args.downloadLibrariesType, asperaKey, sample_outdir, args.downloadCramBam, args.threads, args.downloadInstrumentPlatform)
        else:
            fastq_files = listIDs[sample]

        fileSize = None

        run_successfully_rematch_first = None
        run_successfully_rematch_second = None
        time_taken_rematch_first = 0
        time_taken_rematch_second = 0
        if run_successfully_fastq is not False:
            fileSize = sum(os.path.getsize(fastq) for fastq in fastq_files)
            # Run ReMatCh
            time_taken_rematch_first, run_successfully_rematch_first, data_by_gene, sample_data_general_first, consensus_files, consensus_sequences = rematch_module.runRematchModule(sample, fastq_files, reference_file, args.threads, sample_outdir, args.extraSeq, args.minCovPresence, args.minCovCall, args.minFrequencyDominantAllele, args.minGeneCoverage, args.conservedSeq, args.debug, args.numMapLoc, args.minGeneIdentity, 'first', args.softClip_baseQuality, args.softClip_recodeRun, reference_dict, args.softClip_cigarFlagRecode, args.bowtieOPT, gene_list_reference, args.notWriteConsensus)
            if run_successfully_rematch_first:
                if args.mlst is not None and (args.mlstRun == 'first' or args.mlstRun == 'all'):
                    run_get_st(sample, mlst_dicts, consensus_sequences, args.mlstConsensus, 'first', workdir, time_str)
                write_data_by_gene(gene_list_reference, args.minGeneCoverage, sample, data_by_gene, workdir, time_str, 'first_run', args.minGeneIdentity, 'coverage_depth')
                if args.reportSequenceCoverage:
                    write_data_by_gene(gene_list_reference, args.minGeneCoverage, sample, data_by_gene, workdir, time_str, 'first_run', args.minGeneIdentity, 'sequence_coverage')
                if args.doubleRun:
                    rematch_second_outdir = os.path.join(sample_outdir, 'rematch_second_run', '')
                    if not os.path.isdir(rematch_second_outdir):
                        os.mkdir(rematch_second_outdir)
                    consensus_concatenated_fasta, consensus_concatenated_gene_list, consensus_concatenated_dict = concatenate_extraSeq_2_consensus(consensus_files['noMatter'], reference_file, args.extraSeq, rematch_second_outdir)
                    if len(consensus_concatenated_gene_list) > 0:
                        time_taken_rematch_second, run_successfully_rematch_second, data_by_gene, sample_data_general_second, consensus_files, consensus_sequences = rematch_module.runRematchModule(sample, fastq_files, consensus_concatenated_fasta, args.threads, rematch_second_outdir, args.extraSeq, args.minCovPresence, args.minCovCall, args.minFrequencyDominantAllele, args.minGeneCoverage, args.conservedSeq, args.debug, args.numMapLoc, args.minGeneIdentity, 'second', args.softClip_baseQuality, args.softClip_recodeRun, consensus_concatenated_dict, args.softClip_cigarFlagRecode, args.bowtieOPT, gene_list_reference, args.notWriteConsensus)
                        if not args.debug:
                            os.remove(consensus_concatenated_fasta)
                        if run_successfully_rematch_second:
                            if args.mlst is not None and (args.mlstRun == 'second' or args.mlstRun == 'all'):
                                run_get_st(sample, mlst_dicts, consensus_sequences, args.mlstConsensus, 'second', workdir, time_str)
                            write_data_by_gene(gene_list_reference, args.minGeneCoverage, sample, data_by_gene, workdir, time_str, 'second_run', args.minGeneIdentity, 'coverage_depth')
                            if args.reportSequenceCoverage:
                                write_data_by_gene(gene_list_reference, args.minGeneCoverage, sample, data_by_gene, workdir, time_str, 'second_run', args.minGeneIdentity, 'sequence_coverage')
                    else:
                        print 'No sequences left after ReMatCh module first run. Second run will not be performed'
                        if os.path.isfile(consensus_concatenated_fasta):
                            os.remove(consensus_concatenated_fasta)
                        if os.path.isdir(rematch_second_outdir):
                            utils.removeDirectory(rematch_second_outdir)

        if not searched_fastq_files and not args.keepDownloadedFastq and fastq_files is not None:
            for fastq in fastq_files:
                if os.path.isfile(fastq):
                    os.remove(fastq)

        time_taken = utils.runTime(sample_start_time)

        write_sample_report(sample, workdir, time_str, fileSize, run_successfully_fastq, run_successfully_rematch_first, run_successfully_rematch_second, time_taken_fastq, time_taken_rematch_first, time_taken_rematch_second, time_taken, sequencingInformation, sample_data_general_first if run_successfully_rematch_first else {'number_absent_genes': None, 'number_genes_multiple_alleles': None, 'mean_sample_coverage': None}, sample_data_general_second if run_successfully_rematch_second else {'number_absent_genes': None, 'number_genes_multiple_alleles': None, 'mean_sample_coverage': None}, fastq_files if fastq_files is not None else '')

        if all([run_successfully_fastq is not False, run_successfully_rematch_first is not False, run_successfully_rematch_second is not False]):
            number_samples_successfully += 1

    return number_samples_successfully, len(listIDs)


def main():
    parser = argparse.ArgumentParser(prog='rematch.py', description='Reads mapping against target sequences, checking mapping and consensus sequences production', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-r', '--reference', type=argparse.FileType('r'), metavar='/path/to/reference_sequence.fasta', help='Fasta file containing reference sequences', required=False)
    parser_optional_general.add_argument('-w', '--workdir', type=str, metavar='/path/to/workdir/directory/', help='Path to the directory where ReMatCh will run and produce the outputs with reads (ended with fastq.gz/fq.gz and, in case of PE data, pair-end direction coded as _R1_001 / _R2_001 or _1 / _2) already present (organized in sample folders) or to be downloaded', required=False, default='.')
    parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads to use', required=False, default=1)
    parser_optional_general.add_argument('--mlst', type=str, metavar='"Streptococcus agalactiae"', help='Species name (same as in PubMLST) to be used in MLST determination. ReMatCh will use Bowtie2 very-sensitive-local mapping parameters and will recode the soft clip CIGAR flags of the first run', required=False)
    parser_optional_general.add_argument('--doNotUseProvidedSoftware', action='store_true', help='Tells ReMatCh to not use Bowtie2, Samtools and Bcftools that are provided with it')

    parser_optional_rematch = parser.add_argument_group('ReMatCh module facultative options')
    parser_optional_rematch.add_argument('--conservedSeq', action='store_true', help=argparse.SUPPRESS)
    # parser_optional_rematch.add_argument('--conservedSeq', action='store_true', help='This option can be used with conserved sequences like MLST genes to speedup the analysis by alignning reads using Bowtie2 sensitive algorithm')
    parser_optional_rematch.add_argument('--extraSeq', type=int, metavar='N', help='Sequence length added to both ends of target sequences (usefull to improve reads mapping to the target one) that will be trimmed in ReMatCh outputs', required=False, default=0)
    parser_optional_rematch.add_argument('--minCovPresence', type=int, metavar='N', help='Reference position minimum coverage depth to consider the position to be present in the sample', required=False, default=5)
    parser_optional_rematch.add_argument('--minCovCall', type=int, metavar='N', help='Reference position minimum coverage depth to perform a base call. Lower coverage will be coded as N', required=False, default=10)
    parser_optional_rematch.add_argument('--minFrequencyDominantAllele', type=float, metavar='0.6', help='Minimum relative frequency of the dominant allele coverage depth (value between [0, 1]). Positions with lower values will be considered as having multiple alleles (and will be coded as N)', required=False, default=0.6)
    parser_optional_rematch.add_argument('--minGeneCoverage', type=int, metavar='N', help='Minimum percentage of target reference gene sequence covered by --minCovPresence to consider a gene to be present (value between [0, 100])', required=False, default=70)
    parser_optional_rematch.add_argument('--minGeneIdentity', type=int, metavar='N', help='Minimum percentage of identity of reference gene sequence covered by --minCovCall to consider a gene to be present (value between [0, 100]). One INDEL will be considered as one difference', required=False, default=80)
    parser_optional_rematch.add_argument('--numMapLoc', type=int, metavar='N', help=argparse.SUPPRESS, required=False, default=1)
    # parser_optional_rematch.add_argument('--numMapLoc', type=int, metavar='N', help='Maximum number of locations to which a read can map (sometimes useful when mapping against similar sequences)', required=False, default=1)
    parser_optional_rematch.add_argument('--doubleRun', action='store_true', help='Tells ReMatCh to run a second time using as reference the noMatter consensus sequence produced in the first run. This will improve consensus sequence determination for sequences with high percentage of target reference gene sequence covered')
    parser_optional_rematch.add_argument('--reportSequenceCoverage', action='store_true', help='Produce an extra combined_report.data_by_gene with the sequence coverage instead of coverage depth')
    parser_optional_rematch.add_argument('--notWriteConsensus', action='store_true', help='Do not write consensus sequences')
    parser_optional_rematch.add_argument('--bowtieOPT', type=str, metavar='"--no-mixed"', help='Extra Bowtie2 options', required=False)
    parser_optional_rematch.add_argument('--debug', action='store_true', help='DeBug Mode: do not remove temporary files')

    parser_optional_mlst = parser.add_argument_group('MLST facultative options')
    parser_optional_rematch.add_argument('--mlstReference', action='store_true', help='If the curated scheme for MLST alleles is available, tells ReMatCh to use these as reference (force Bowtie2 to run with very-sensitive-local parameters, and sets --extraSeq to 200), otherwise ReMatCh uses the first alleles of each MLST gene fragment in PubMLST as reference sequences (force Bowtie2 to run with very-sensitive-local parameters, and sets --extraSeq to 0)')
    parser_optional_mlst.add_argument('--mlstSchemaNumber', type=int, metavar='N', help='Number of the species PubMLST schema to be used in case of multiple schemes available (by default will use the first schema)', required=False)
    parser_optional_mlst.add_argument('--mlstConsensus', choices=['noMatter', 'correct', 'alignment', 'all'], type=str, metavar='noMatter', help='Consensus sequence to be used in MLST determination (available options: %(choices)s)', required=False, default='noMatter')
    parser_optional_mlst.add_argument('--mlstRun', choices=['first', 'second', 'all'], type=str, metavar='first', help='ReMatCh run outputs to be used in MLST determination (available options: %(choices)s)', required=False, default='all')

    parser_optional_download = parser.add_argument_group('Download facultative options')
    parser_optional_download.add_argument('-a', '--asperaKey', type=argparse.FileType('r'), metavar='/path/to/asperaweb_id_dsa.openssh', help='Tells ReMatCh to download fastq files from ENA using Aspera Connect. With this option, the path to Private-key file asperaweb_id_dsa.openssh must be provided (normaly found in ~/.aspera/connect/etc/asperaweb_id_dsa.openssh).', required=False)
    parser_optional_download.add_argument('-k', '--keepDownloadedFastq', action='store_true', help='Tells ReMatCh to keep the fastq files downloaded')
    parser_optional_download.add_argument('--downloadLibrariesType', type=str, metavar='PAIRED', help='Tells ReMatCh to download files with specific library layout (available options: %(choices)s)', choices=['PAIRED', 'SINGLE', 'BOTH'], required=False, default='BOTH')
    parser_optional_download.add_argument('--downloadInstrumentPlatform', type=str, metavar='ILLUMINA', help='Tells ReMatCh to download files with specific library layout (available options: %(choices)s)', choices=['ILLUMINA', 'ALL'], required=False, default='ILLUMINA')
    parser_optional_download.add_argument('--downloadCramBam', action='store_true', help='Tells ReMatCh to also download cram/bam files and convert them to fastq files')

    parser_optional_softClip = parser.add_argument_group('Soft clip facultative options')
    parser_optional_softClip.add_argument('--softClip_baseQuality', type=int, metavar='N', help='Base quality phred score in reads soft clipped regions', required=False, default=7)
    parser_optional_download.add_argument('--softClip_recodeRun', type=str, metavar='first', help='ReMatCh run to recode soft clipped regions (available options: %(choices)s)', choices=['first', 'second', 'both', 'none'], required=False, default='none')
    parser_optional_download.add_argument('--softClip_cigarFlagRecode', type=str, metavar='M', help='CIGAR flag to recode CIGAR soft clip (available options: %(choices)s)', choices=['M', 'I', 'X'], required=False, default='X')

    parser_optional_download_exclusive = parser.add_mutually_exclusive_group()
    parser_optional_download_exclusive.add_argument('-l', '--listIDs', type=argparse.FileType('r'), metavar='/path/to/list_IDs.txt', help='Path to list containing the IDs to be downloaded (one per line)', required=False)
    parser_optional_download_exclusive.add_argument('-t', '--taxon', type=str, metavar='"Streptococcus agalactiae"', help='Taxon name for which ReMatCh will download fastq files', required=False)

    args = parser.parse_args()

    if args.reference is None and not args.mlstReference:
        parser.error('At least --reference or --mlstReference should be provided')
    elif args.reference is not None and args.mlstReference:
        parser.error('Only --reference or --mlstReference should be provided')
    else:
        if args.mlstReference:
            if args.mlst is None:
                parser.error('Please provide species name using --mlst')

    if args.minFrequencyDominantAllele < 0 or args.minFrequencyDominantAllele > 1:
        parser.error('--minFrequencyDominantAllele should be a value between [0, 1]')

    if args.minGeneCoverage < 0 or args.minGeneCoverage > 100:
        parser.error('--minGeneCoverage should be a value between [0, 100]')
    if args.minGeneIdentity < 0 or args.minGeneIdentity > 100:
        parser.error('--minGeneIdentity should be a value between [0, 100]')

    start_time = time.time()

    number_samples_successfully, samples_total_number = runRematch(args)

    print '\n' + 'END ReMatCh'
    print '\n' + str(number_samples_successfully) + ' samples out of ' + str(samples_total_number) + ' run successfully'
    time_taken = utils.runTime(start_time)
    del time_taken

    if number_samples_successfully == 0:
        sys.exit('No samples run successfully!')


if __name__ == "__main__":
    main()

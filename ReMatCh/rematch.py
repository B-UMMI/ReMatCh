#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
rematch.py - Reads mapping against target sequences, checking mapping
and consensus sequences production
<https://github.com/B-UMMI/ReMatCh/>

Copyright (C) 2019 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: August 08, 2019

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

try:
    from __init__ import __version__

    import modules.utils as utils
    import modules.seqFromWebTaxon as seq_from_web_taxon
    import modules.download as download
    import modules.rematch_module as rematch_module
    import modules.checkMLST as check_mlst
except ImportError:
    from ReMatCh.__init__ import __version__

    from ReMatCh.modules import utils as utils
    from ReMatCh.modules import seqFromWebTaxon as seq_from_web_taxon
    from ReMatCh.modules import download as download
    from ReMatCh.modules import rematch_module as rematch_module
    from ReMatCh.modules import checkMLST as check_mlst


def search_fastq_files(directory):
    files_extensions = ['.fastq.gz', '.fq.gz']
    pair_end_files_separation = [['_R1_001.f', '_R2_001.f'], ['_1.f', '_2.f']]

    list_ids = {}
    directories = [d for d in os.listdir(directory) if
                   not d.startswith('.') and os.path.isdir(os.path.join(directory, d, ''))]
    for directory_found in directories:
        if directory_found != 'pubmlst':
            directory_path = os.path.join(directory, directory_found, '')

            fastq_found = []
            files = [f for f in os.listdir(directory_path) if
                     not f.startswith('.') and os.path.isfile(os.path.join(directory_path, f))]
            for file_found in files:
                if file_found.endswith(tuple(files_extensions)):
                    fastq_found.append(file_found)

            if len(fastq_found) == 1:
                list_ids[directory_found] = [os.path.join(directory_path, f) for f in fastq_found]
            elif len(fastq_found) >= 2:
                file_pair = []

                # Search pairs
                for pe_separation in pair_end_files_separation:
                    for fastq in fastq_found:
                        if pe_separation[0] in fastq or pe_separation[1] in fastq:
                            file_pair.append(fastq)

                    if len(file_pair) == 2:
                        break
                    else:
                        file_pair = []

                # Search single
                if len(file_pair) == 0:
                    for pe_separation in pair_end_files_separation:
                        for fastq in fastq_found:
                            if pe_separation[0] not in fastq or pe_separation[1] not in fastq:
                                file_pair.append(fastq)

                    if len(file_pair) >= 1:
                        file_pair = file_pair[0]

                if len(file_pair) >= 1:
                    list_ids[directory_found] = [os.path.join(directory_path, f) for f in file_pair]

    return list_ids


def get_list_ids_from_file(file_list_ids):
    list_ids = []

    with open(file_list_ids, 'rtU') as lines:
        for line in lines:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                list_ids.append(line)

    if len(list_ids) == 0:
        sys.exit('No runIDs were found in ' + file_list_ids)

    return list_ids


def get_taxon_run_ids(taxon_name, outputfile):
    seq_from_web_taxon.run_seq_from_web_taxon(taxon_name, outputfile, True, True, True, False)

    run_ids = []
    with open(outputfile, 'rtU') as reader:
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                if not line.startswith('#'):
                    line = line.split('\t')
                    run_ids.append(line[0])

    return run_ids


def get_list_ids(workdir, file_list_ids, taxon_name):
    searched_fastq_files = False
    list_ids = []
    if file_list_ids is None and taxon_name is None:
        list_ids = search_fastq_files(workdir)
        searched_fastq_files = True
    elif file_list_ids is not None:
        list_ids = get_list_ids_from_file(os.path.abspath(file_list_ids))
    elif taxon_name is not None and file_list_ids is None:
        list_ids = get_taxon_run_ids(taxon_name, os.path.join(workdir, 'IDs_list.seqFromWebTaxon.tab'))

    if len(list_ids) == 0:
        sys.exit('No IDs were found')
    return list_ids, searched_fastq_files


def format_gene_info(gene_specific_info, minimum_gene_coverage, minimum_gene_identity, reported_data_type, summary,
                     sample, genes_present):
    info = None
    if gene_specific_info['gene_coverage'] >= minimum_gene_coverage and \
            gene_specific_info['gene_identity'] >= minimum_gene_identity:
        if summary and sample not in genes_present:
            genes_present[sample] = {}

        if gene_specific_info['gene_number_positions_multiple_alleles'] == 0:
            s = str(gene_specific_info[reported_data_type])
            info = str(s)
            if summary:
                genes_present[sample][gene_specific_info['header']] = str(s)
        else:
            s = 'multiAlleles_' + str(gene_specific_info[reported_data_type])
            info = str(s)
            if summary:
                genes_present[sample][gene_specific_info['header']] = str(s)
    else:
        info = 'absent_' + str(gene_specific_info[reported_data_type])

    return info, genes_present


def write_data_by_gene(gene_list_reference, minimum_gene_coverage, sample, data_by_gene, outdir, time_str, run_times,
                       minimum_gene_identity, reported_data_type, summary, genes_present):
    combined_report = \
        os.path.join(outdir,
                     'combined_report.data_by_gene.' + run_times + '.' + reported_data_type + '.' + time_str + '.tab')

    if reported_data_type == 'coverage_depth':
        reported_data_type = 'gene_mean_read_coverage'
    elif reported_data_type == 'sequence_coverage':
        reported_data_type = 'gene_coverage'

    combined_report_exist = os.path.isfile(combined_report)
    with open(combined_report, 'at') as writer:
        seq_list = sorted(gene_list_reference.keys())
        if not combined_report_exist:
            writer.write('#sample' + '\t' + '\t'.join([gene_list_reference[seq] for seq in seq_list]) + '\n')

        results = {}
        headers = []
        for i in data_by_gene:
            results[data_by_gene[i]['header']], genes_present = format_gene_info(data_by_gene[i], minimum_gene_coverage,
                                                                                 minimum_gene_identity,
                                                                                 reported_data_type, summary, sample,
                                                                                 genes_present)
            headers.append(data_by_gene[i]['header'])

        if len(headers) != gene_list_reference:
            for gene in gene_list_reference:
                if gene not in headers:
                    results[gene] = 'NA'

        writer.write(sample + '\t' + '\t'.join([results[seq] for seq in seq_list]) + '\n')

    return genes_present


def write_sample_report(sample, outdir, time_str, file_size, run_successfully_fastq, run_successfully_rematch_first,
                        run_successfully_rematch_second, time_taken_fastq, time_taken_rematch_first,
                        time_taken_rematch_second, time_taken_sample, sequencing_information, sample_data_general_first,
                        sample_data_general_second, fastq_used):
    sample_report = os.path.join(outdir, 'sample_report.' + time_str + '.tab')
    report_exist = os.path.isfile(sample_report)

    header_general = ['sample', 'sample_run_successfully', 'sample_run_time', 'files_size', 'download_run_successfully',
                      'download_run_time', 'rematch_run_successfully_first', 'rematch_run_time_first',
                      'rematch_run_successfully_second', 'rematch_run_time_second']
    header_data_general = ['number_absent_genes', 'number_genes_multiple_alleles', 'mean_sample_coverage']
    header_sequencing = ['run_accession', 'instrument_platform', 'instrument_model', 'library_layout', 'library_source',
                         'extra_run_accession', 'nominal_length', 'read_count', 'base_count', 'date_download']

    with open(sample_report, 'at') as writer:
        if not report_exist:
            writer.write('#' + '\t'.join(header_general) + '\t' + '_first\t'.join(header_data_general) + '_first\t' +
                         '_second\t'.join(header_data_general) + '_second\t' + '\t'.join(header_sequencing) + '\t' +
                         'fastq_used' + '\n')

        writer.write('\t'.join([sample,
                                str(all([run_successfully_fastq is not False,
                                         run_successfully_rematch_first is not False,
                                         run_successfully_rematch_second is not False])),
                                str(time_taken_sample),
                                str(file_size),
                                str(run_successfully_fastq),
                                str(time_taken_fastq),
                                str(run_successfully_rematch_first),
                                str(time_taken_rematch_first),
                                str(run_successfully_rematch_second),
                                str(time_taken_rematch_second)]) +
                     '\t' + '\t'.join([str(sample_data_general_first[i]) for i in header_data_general]) +
                     '\t' + '\t'.join([str(sample_data_general_second[i]) for i in header_data_general]) +
                     '\t' + '\t'.join([str(sequencing_information[i]) for i in header_sequencing]) +
                     '\t' + ','.join(fastq_used) + '\n')


def concatenate_extra_seq_2_consensus(consensus_sequence, reference_sequence, extra_seq_length, outdir):
    reference_dict, ignore, ignore = rematch_module.get_sequence_information(reference_sequence, extra_seq_length)
    consensus_dict, genes, ignore = rematch_module.get_sequence_information(consensus_sequence, 0)
    number_consensus_with_sequences = 0
    for k, values_consensus in list(consensus_dict.items()):
        for values_reference in list(reference_dict.values()):
            if values_reference['header'] == values_consensus['header']:
                if len(set(consensus_dict[k]['sequence'])) > 1:
                    number_consensus_with_sequences += 1
                    if extra_seq_length <= len(values_reference['sequence']):
                        right_extra_seq = \
                            '' if extra_seq_length == 0 else values_reference['sequence'][-extra_seq_length:]
                        consensus_dict[k]['sequence'] = \
                            values_reference['sequence'][:extra_seq_length] + \
                            consensus_dict[k]['sequence'] + \
                            right_extra_seq
                        consensus_dict[k]['length'] += extra_seq_length + len(right_extra_seq)

    consensus_concatenated = os.path.join(outdir, 'consensus_concatenated_extraSeq.fasta')
    with open(consensus_concatenated, 'wt') as writer:
        for i in consensus_dict:
            writer.write('>' + consensus_dict[i]['header'] + '\n')
            fasta_sequence_lines = rematch_module.chunkstring(consensus_dict[i]['sequence'], 80)
            for line in fasta_sequence_lines:
                writer.write(line + '\n')

    return consensus_concatenated, genes, consensus_dict, number_consensus_with_sequences


def clean_headers_reference_file(reference_file, outdir, extra_seq):
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", ":"]
    print('Checking if reference sequences contain ' + str(problematic_characters) + '\n')
    # headers_changed = False
    new_reference_file = str(reference_file)
    sequences, genes, headers_changed = rematch_module.get_sequence_information(reference_file, extra_seq)
    if headers_changed:
        print('At least one of the those characters was found. Replacing those with _' + '\n')
        new_reference_file = \
            os.path.join(outdir, os.path.splitext(os.path.basename(reference_file))[0] + '.headers_renamed.fasta')
        with open(new_reference_file, 'wt') as writer:
            for i in sequences:
                writer.write('>' + sequences[i]['header'] + '\n')
                fasta_sequence_lines = rematch_module.chunkstring(sequences[i]['sequence'], 80)
                for line in fasta_sequence_lines:
                    writer.write(line + '\n')
    return new_reference_file, genes, sequences


def write_mlst_report(sample, run_times, consensus_type, st, alleles_profile, loci_order, outdir, time_str):
    mlst_report = os.path.join(outdir, 'mlst_report.' + time_str + '.tab')
    mlst_report_exist = os.path.isfile(mlst_report)
    with open(mlst_report, 'at') as writer:
        if not mlst_report_exist:
            writer.write('\t'.join(['#sample', 'ReMatCh_run', 'consensus_type', 'ST'] + loci_order) + '\n')
        writer.write('\t'.join([sample, run_times, consensus_type, str(st)] + alleles_profile.split(',')) + '\n')


def run_get_st(sample, mlst_dicts, consensus_sequences, mlst_consensus, run_times, outdir, time_str):
    if mlst_consensus == 'all':
        for consensus_type in consensus_sequences:
            print('Searching MLST for ' + consensus_type + ' consensus')
            st, alleles_profile = check_mlst.get_st(mlst_dicts, consensus_sequences[consensus_type])
            write_mlst_report(sample, run_times, consensus_type, st, alleles_profile, mlst_dicts[2], outdir, time_str)
            print('ST found: ' + str(st) + ' (' + alleles_profile + ')')
    else:
        st, alleles_profile = check_mlst.get_st(mlst_dicts, consensus_sequences[mlst_consensus])
        write_mlst_report(sample, run_times, mlst_consensus, st, alleles_profile, mlst_dicts[2], outdir, time_str)
        print('ST found for ' + mlst_consensus + ' consensus: ' + str(st) + ' (' + alleles_profile + ')')


def write_summary_report(outdir, reported_data_type, time_str, gene_list_reference, genes_present):
    with open(os.path.join(outdir,
                           'summary.{reported_data_type}.{time_str}.tab'.format(reported_data_type=reported_data_type,
                                                                                time_str=time_str)), 'wt') as writer:
        seq_list = []
        for info in list(genes_present.values()):
            seq_list.extend(list(info.keys()))
            seq_list = list(set(seq_list))
        writer.write('#sample' + '\t' + '\t'.join([gene_list_reference[seq] for seq in sorted(seq_list)]) + '\n')
        for sample, info in list(genes_present.items()):
            data = []
            for seq in sorted(seq_list):
                if seq in info:
                    data.append(info[seq])
                else:
                    data.append('NF')
            writer.write(sample + '\t' + '\t'.join(data) + '\n')


def run_rematch(args):
    workdir = os.path.abspath(args.workdir)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    aspera_key = os.path.abspath(args.asperaKey.name) if args.asperaKey is not None else None

    # Start logger
    logfile, time_str = utils.start_logger(workdir)

    # Get general information
    script_path = utils.general_information(logfile, __version__, workdir, time_str, args.doNotUseProvidedSoftware,
                                            aspera_key, args.downloadCramBam, args.SRA, args.SRAopt)

    # Set list_ids
    list_ids, searched_fastq_files = get_list_ids(workdir, args.listIDs.name if args.listIDs is not None else None,
                                                  args.taxon)

    mlst_sequences = None
    mlst_dicts = None
    if args.mlst is not None:
        time_taken_pub_mlst, mlst_dicts, mlst_sequences = check_mlst.download_pub_mlst_xml(args.mlst,
                                                                                           args.mlstSchemaNumber,
                                                                                           workdir)
        args.softClip_recodeRun = 'first'

    if args.reference is None:
        if args.mlst is not None:
            reference_file = check_mlst.check_existing_schema(args.mlst, args.mlstSchemaNumber, script_path)
            args.extraSeq = 200
            if reference_file is None:
                print('It was not found provided MLST scheme sequences for ' + args.mlst)
                print('Trying to obtain reference MLST sequences from PubMLST')
                if len(mlst_sequences) > 0:
                    reference_file = check_mlst.write_mlst_reference(args.mlst, mlst_sequences, workdir, time_str)
                    args.extraSeq = 0
                else:
                    sys.exit('It was not possible to download MLST sequences from PubMLST!')
            else:
                print('Using provided scheme as referece: ' + reference_file)
        else:
            sys.exit('Need to provide at least one of the following options: "--reference" and "--mlst"')
    else:
        reference_file = os.path.abspath(args.reference.name)

    # Run ReMatCh for each sample
    print('\n' + 'STARTING ReMatCh' + '\n')

    # Clean sequences headers
    reference_file, gene_list_reference, reference_dict = clean_headers_reference_file(reference_file, workdir,
                                                                                       args.extraSeq)

    if args.mlst is not None:
        problem_genes = False
        for header in mlst_sequences:
            if header not in gene_list_reference:
                print('MLST gene {header} not found between reference sequences'.format(header=header))
                problem_genes = True
        if problem_genes:
            sys.exit('Missing MLST genes from reference sequences (at least sequences names do not match)!')

    if len(gene_list_reference) == 0:
        sys.exit('No sequences left')

    # To use in combined report

    number_samples_successfully = 0
    genes_present_coverage_depth = {}
    genes_present_sequence_coverage = {}
    for sample in list_ids:
        sample_start_time = time.time()
        print('\n\n' + 'Sample ID: ' + sample)

        # Create sample outdir
        sample_outdir = os.path.join(workdir, sample, '')
        if not os.path.isdir(sample_outdir):
            os.mkdir(sample_outdir)

        run_successfully_fastq = None
        time_taken_fastq = 0
        sequencing_information = {'run_accession': None, 'instrument_platform': None, 'instrument_model': None,
                                  'library_layout': None, 'library_source': None, 'extra_run_accession': None,
                                  'nominal_length': None, 'read_count': None, 'base_count': None, 'date_download': None}
        if not searched_fastq_files:
            # Download Files
            time_taken_fastq, run_successfully_fastq, fastq_files, sequencing_information = \
                download.run_download(sample, args.downloadLibrariesType, aspera_key, sample_outdir,
                                      args.downloadCramBam, args.threads, args.downloadInstrumentPlatform, args.SRA,
                                      args.SRAopt)
        else:
            fastq_files = list_ids[sample]

        file_size = None

        run_successfully_rematch_first = None
        run_successfully_rematch_second = None
        time_taken_rematch_first = 0
        time_taken_rematch_second = 0
        sample_data_general_first = None
        sample_data_general_second = None
        if run_successfully_fastq is not False:
            file_size = sum(os.path.getsize(fastq) for fastq in fastq_files)
            # Run ReMatCh
            time_taken_rematch_first, run_successfully_rematch_first, data_by_gene, sample_data_general_first, \
            consensus_files, consensus_sequences = \
                rematch_module.run_rematch_module(sample, fastq_files, reference_file, args.threads, sample_outdir,
                                                  args.extraSeq, args.minCovPresence, args.minCovCall,
                                                  args.minFrequencyDominantAllele, args.minGeneCoverage,
                                                  args.debug, args.numMapLoc, args.minGeneIdentity,
                                                  'first', args.softClip_baseQuality, args.softClip_recodeRun,
                                                  reference_dict, args.softClip_cigarFlagRecode,
                                                  args.bowtieAlgo, args.bowtieOPT,
                                                  gene_list_reference, args.notWriteConsensus, clean_run=True)
            if run_successfully_rematch_first:
                if args.mlst is not None and (args.mlstRun == 'first' or args.mlstRun == 'all'):
                    run_get_st(sample, mlst_dicts, consensus_sequences, args.mlstConsensus, 'first', workdir, time_str)
                genes_present_coverage_depth = write_data_by_gene(gene_list_reference, args.minGeneCoverage, sample,
                                                                  data_by_gene, workdir, time_str, 'first_run',
                                                                  args.minGeneIdentity, 'coverage_depth', args.summary,
                                                                  genes_present_coverage_depth)
                if args.reportSequenceCoverage:
                    genes_present_sequence_coverage = write_data_by_gene(gene_list_reference, args.minGeneCoverage,
                                                                         sample, data_by_gene, workdir, time_str,
                                                                         'first_run', args.minGeneIdentity,
                                                                         'sequence_coverage', args.summary,
                                                                         genes_present_sequence_coverage)
                if args.doubleRun:
                    rematch_second_outdir = os.path.join(sample_outdir, 'rematch_second_run', '')
                    if not os.path.isdir(rematch_second_outdir):
                        os.mkdir(rematch_second_outdir)
                    consensus_concatenated_fasta, consensus_concatenated_gene_list, consensus_concatenated_dict, \
                    number_consensus_with_sequences = \
                        concatenate_extra_seq_2_consensus(consensus_files['noMatter'], reference_file, args.extraSeq,
                                                          rematch_second_outdir)
                    if len(consensus_concatenated_gene_list) > 0:
                        if args.mlst is None or \
                                (args.mlst is not None and number_consensus_with_sequences == len(gene_list_reference)):
                            time_taken_rematch_second, run_successfully_rematch_second, data_by_gene, \
                            sample_data_general_second, consensus_files, consensus_sequences = \
                                rematch_module.run_rematch_module(sample, fastq_files, consensus_concatenated_fasta,
                                                                  args.threads, rematch_second_outdir, args.extraSeq,
                                                                  args.minCovPresence, args.minCovCall,
                                                                  args.minFrequencyDominantAllele, args.minGeneCoverage,
                                                                  args.debug, args.numMapLoc,
                                                                  args.minGeneIdentity, 'second',
                                                                  args.softClip_baseQuality, args.softClip_recodeRun,
                                                                  consensus_concatenated_dict,
                                                                  args.softClip_cigarFlagRecode,
                                                                  args.bowtieAlgo, args.bowtieOPT,
                                                                  gene_list_reference, args.notWriteConsensus,
                                                                  clean_run=True)
                            if not args.debug:
                                os.remove(consensus_concatenated_fasta)
                            if run_successfully_rematch_second:
                                if args.mlst is not None and (args.mlstRun == 'second' or args.mlstRun == 'all'):
                                    run_get_st(sample, mlst_dicts, consensus_sequences, args.mlstConsensus, 'second',
                                               workdir, time_str)
                                _ = write_data_by_gene(gene_list_reference, args.minGeneCoverage, sample, data_by_gene,
                                                       workdir, time_str, 'second_run', args.minGeneIdentity,
                                                       'coverage_depth', False, {})
                                if args.reportSequenceCoverage:
                                    _ = write_data_by_gene(gene_list_reference, args.minGeneCoverage, sample,
                                                           data_by_gene, workdir, time_str, 'second_run',
                                                           args.minGeneIdentity, 'sequence_coverage', False, {})
                        else:
                            print('Some sequences missing after ReMatCh module first run. Second run will not be'
                                  ' performed')
                            if os.path.isfile(consensus_concatenated_fasta):
                                os.remove(consensus_concatenated_fasta)
                            if os.path.isdir(rematch_second_outdir):
                                utils.remove_directory(rematch_second_outdir)
                    else:
                        print('No sequences left after ReMatCh module first run. Second run will not be performed')
                        if os.path.isfile(consensus_concatenated_fasta):
                            os.remove(consensus_concatenated_fasta)
                        if os.path.isdir(rematch_second_outdir):
                            utils.remove_directory(rematch_second_outdir)

        if not searched_fastq_files and not args.keepDownloadedFastq and fastq_files is not None:
            for fastq in fastq_files:
                if os.path.isfile(fastq):
                    os.remove(fastq)

        time_taken = utils.run_time(sample_start_time)

        write_sample_report(sample, workdir, time_str, file_size, run_successfully_fastq,
                            run_successfully_rematch_first, run_successfully_rematch_second, time_taken_fastq,
                            time_taken_rematch_first, time_taken_rematch_second, time_taken, sequencing_information,
                            sample_data_general_first if run_successfully_rematch_first else
                            {'number_absent_genes': None, 'number_genes_multiple_alleles': None,
                             'mean_sample_coverage': None},
                            sample_data_general_second if run_successfully_rematch_second else
                            {'number_absent_genes': None, 'number_genes_multiple_alleles': None,
                             'mean_sample_coverage': None},
                            fastq_files if fastq_files is not None else '')

        if all([run_successfully_fastq is not False,
                run_successfully_rematch_first is not False,
                run_successfully_rematch_second is not False]):
            number_samples_successfully += 1

    if args.summary:
        write_summary_report(workdir, 'coverage_depth', time_str, gene_list_reference, genes_present_coverage_depth)
        if args.reportSequenceCoverage:
            write_summary_report(workdir, 'sequence_coverage', time_str, gene_list_reference,
                                 genes_present_sequence_coverage)

    return number_samples_successfully, len(list_ids)


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 rematch.py"')

    parser = argparse.ArgumentParser(prog='rematch.py',
                                     description='Reads mapping against target sequences, checking mapping and'
                                                 ' consensus sequences production',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version',
                        version='{prog} v{version}'.format(prog=parser.prog, version=__version__))

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-r', '--reference', type=argparse.FileType('r'),
                                         metavar='/path/to/reference_sequence.fasta',
                                         help='Fasta file containing reference sequences', required=False)
    parser_optional_general.add_argument('-w', '--workdir', type=str, metavar='/path/to/workdir/directory/',
                                         help='Path to the directory where ReMatCh will run and produce the outputs'
                                              ' with reads (ended with fastq.gz/fq.gz and, in case of PE data, pair-end'
                                              ' direction coded as _R1_001 / _R2_001 or _1 / _2) already'
                                              ' present (organized in sample folders) or to be downloaded',
                                         required=False, default='.')
    parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads to use',
                                         required=False, default=1)
    parser_optional_general.add_argument('--mlst', type=str, metavar='"Streptococcus agalactiae"',
                                         help='Species name (same as in PubMLST) to be used in MLST'
                                              ' determination. ReMatCh will use Bowtie2 very-sensitive-local mapping'
                                              ' parameters and will recode the soft clip CIGAR flags of the first run',
                                         required=False)
    parser_optional_general.add_argument('--doNotUseProvidedSoftware', action='store_true',
                                         help='Tells ReMatCh to not use Bowtie2, Samtools and Bcftools that are'
                                              ' provided with it')

    parser_optional_download_exclusive = parser.add_mutually_exclusive_group()
    parser_optional_download_exclusive.add_argument('-l', '--listIDs', type=argparse.FileType('r'),
                                                    metavar='/path/to/list_IDs.txt',
                                                    help='Path to list containing the IDs to be'
                                                         ' downloaded (one per line)', required=False)
    parser_optional_download_exclusive.add_argument('-t', '--taxon', type=str, metavar='"Streptococcus agalactiae"',
                                                    help='Taxon name for which ReMatCh will download fastq files',
                                                    required=False)

    parser_optional_rematch = parser.add_argument_group('ReMatCh module facultative options')
    parser_optional_rematch.add_argument('--extraSeq', type=int, metavar='N',
                                         help='Sequence length added to both ends of target sequences (usefull to'
                                              ' improve reads mapping to the target one) that will be trimmed in'
                                              ' ReMatCh outputs', required=False, default=0)
    parser_optional_rematch.add_argument('--minCovPresence', type=int, metavar='N',
                                         help='Reference position minimum coverage depth to consider the position to be'
                                              ' present in the sample', required=False, default=5)
    parser_optional_rematch.add_argument('--minCovCall', type=int, metavar='N',
                                         help='Reference position minimum coverage depth to perform a base call. Lower'
                                              ' coverage will be coded as N', required=False, default=10)
    parser_optional_rematch.add_argument('--minFrequencyDominantAllele', type=float, metavar='0.6',
                                         help='Minimum relative frequency of the dominant allele coverage depth (value'
                                              ' between [0, 1]). Positions with lower values will be considered as'
                                              ' having multiple alleles (and will be coded as N)', required=False,
                                         default=0.6)
    parser_optional_rematch.add_argument('--minGeneCoverage', type=int, metavar='N',
                                         help='Minimum percentage of target reference gene sequence covered'
                                              ' by --minCovPresence to consider a gene to be present (value'
                                              ' between [0, 100])', required=False, default=70)
    parser_optional_rematch.add_argument('--minGeneIdentity', type=int, metavar='N',
                                         help='Minimum percentage of identity of reference gene sequence covered'
                                              ' by --minCovCall to consider a gene to be present (value'
                                              ' between [0, 100]). One INDEL will be considered as one difference',
                                         required=False, default=80)
    parser_optional_rematch.add_argument('--numMapLoc', type=int, metavar='N', help=argparse.SUPPRESS, required=False,
                                         default=1)
    # parser_optional_rematch.add_argument('--numMapLoc', type=int, metavar='N', help='Maximum number of locations to which a read can map (sometimes useful when mapping against similar sequences)', required=False, default=1)
    parser_optional_rematch.add_argument('--doubleRun', action='store_true',
                                         help='Tells ReMatCh to run a second time using as reference the noMatter'
                                              ' consensus sequence produced in the first run. This will improve'
                                              ' consensus sequence determination for sequences with high percentage of'
                                              ' target reference gene sequence covered')
    parser_optional_rematch.add_argument('--reportSequenceCoverage', action='store_true',
                                         help='Produce an extra combined_report.data_by_gene with the sequence coverage'
                                              ' instead of coverage depth')
    parser_optional_rematch.add_argument('--summary', action='store_true',
                                         help='Produce extra report files containing only sequences present in at least'
                                              ' one sample (usefull when using a large number of reference'
                                              ' sequences, and only for first run)')
    parser_optional_rematch.add_argument('--notWriteConsensus', action='store_true',
                                         help='Do not write consensus sequences')
    parser_optional_rematch.add_argument('--bowtieAlgo', type=str, metavar='"--very-sensitive-local"',
                                         help='Bowtie2 alignment mode. It can be an end-to-end alignment (unclipped'
                                              ' alignment) or local alignment (soft clipped alignment). Also, can'
                                              ' choose between fast or sensitive alignments. Please check Bowtie2'
                                              ' manual for extra'
                                              ' information: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml .'
                                              ' This option should be provided between quotes and starting with'
                                              ' an empty space (like --bowtieAlgo " --very-fast") or using equal'
                                              ' sign (like --bowtieAlgo="--very-fast")',
                                         required=False, default='--very-sensitive-local')
    parser_optional_rematch.add_argument('--bowtieOPT', type=str, metavar='"--no-mixed"',
                                         help='Extra Bowtie2 options. This option should be provided between quotes and'
                                              ' starting with an empty space (like --bowtieOPT " --no-mixed") or using'
                                              ' equal sign (like --bowtieOPT="--no-mixed")',
                                         required=False)
    parser_optional_rematch.add_argument('--debug', action='store_true',
                                         help='DeBug Mode: do not remove temporary files')

    parser_optional_mlst = parser.add_argument_group('MLST facultative options')
    parser_optional_rematch.add_argument('--mlstReference', action='store_true',
                                         help='If the curated scheme for MLST alleles is available, tells ReMatCh to'
                                              ' use these as reference (force Bowtie2 to run with very-sensitive-local'
                                              ' parameters, and sets --extraSeq to 200), otherwise ReMatCh uses the'
                                              ' first alleles of each MLST gene fragment in PubMLST as reference'
                                              ' sequences (force Bowtie2 to run with very-sensitive-local parameters,'
                                              ' and sets --extraSeq to 0)')
    parser_optional_mlst.add_argument('--mlstSchemaNumber', type=int, metavar='N',
                                      help='Number of the species PubMLST schema to be used in case of multiple schemes'
                                           ' available (by default will use the first schema)', required=False)
    parser_optional_mlst.add_argument('--mlstConsensus', choices=['noMatter', 'correct', 'alignment', 'all'], type=str,
                                      metavar='noMatter',
                                      help='Consensus sequence to be used in MLST'
                                           ' determination (available options: %(choices)s)', required=False,
                                      default='noMatter')
    parser_optional_mlst.add_argument('--mlstRun', choices=['first', 'second', 'all'], type=str, metavar='first',
                                      help='ReMatCh run outputs to be used in MLST determination (available'
                                           ' options: %(choices)s)', required=False, default='all')

    parser_optional_download = parser.add_argument_group('Download facultative options')
    parser_optional_download.add_argument('-a', '--asperaKey', type=argparse.FileType('r'),
                                          metavar='/path/to/asperaweb_id_dsa.openssh',
                                          help='Tells ReMatCh to download fastq files from ENA using Aspera'
                                               ' Connect. With this option, the path to Private-key file'
                                               ' asperaweb_id_dsa.openssh must be provided (normaly found in'
                                               ' ~/.aspera/connect/etc/asperaweb_id_dsa.openssh).', required=False)
    parser_optional_download.add_argument('-k', '--keepDownloadedFastq', action='store_true',
                                          help='Tells ReMatCh to keep the fastq files downloaded')
    parser_optional_download.add_argument('--downloadLibrariesType', type=str, metavar='PAIRED',
                                          help='Tells ReMatCh to download files with specific library'
                                               ' layout (available options: %(choices)s)',
                                          choices=['PAIRED', 'SINGLE', 'BOTH'], required=False, default='BOTH')
    parser_optional_download.add_argument('--downloadInstrumentPlatform', type=str, metavar='ILLUMINA',
                                          help='Tells ReMatCh to download files with specific library layout (available'
                                               ' options: %(choices)s)', choices=['ILLUMINA', 'ALL'], required=False,
                                          default='ILLUMINA')
    parser_optional_download.add_argument('--downloadCramBam', action='store_true',
                                          help='Tells ReMatCh to also download cram/bam files and convert them to fastq'
                                               ' files')

    parser_optional_sra = parser.add_mutually_exclusive_group()
    parser_optional_sra.add_argument('--SRA', action='store_true',
                                     help='Tells getSeqENA.py to download reads in fastq format only from NCBI SRA'
                                          ' database (not recommended)')
    parser_optional_sra.add_argument('--SRAopt', action='store_true',
                                     help='Tells getSeqENA.py to download reads from NCBI SRA if the download from ENA'
                                          ' fails')

    parser_optional_soft_clip = parser.add_argument_group('Soft clip facultative options')
    parser_optional_soft_clip.add_argument('--softClip_baseQuality', type=int, metavar='N',
                                           help='Base quality phred score in reads soft clipped regions',
                                           required=False,
                                           default=7)
    parser_optional_soft_clip.add_argument('--softClip_recodeRun', type=str, metavar='first',
                                           help='ReMatCh run to recode soft clipped regions (available'
                                                ' options: %(choices)s)', choices=['first', 'second', 'both', 'none'],
                                           required=False, default='none')
    parser_optional_soft_clip.add_argument('--softClip_cigarFlagRecode', type=str, metavar='M',
                                           help='CIGAR flag to recode CIGAR soft clip (available options: %(choices)s)',
                                           choices=['M', 'I', 'X'], required=False, default='X')

    args = parser.parse_args()

    msg = []
    if args.reference is None and not args.mlstReference:
        msg.append('At least --reference or --mlstReference should be provided')
    elif args.reference is not None and args.mlstReference:
        msg.append('Only --reference or --mlstReference should be provided')
    else:
        if args.mlstReference:
            if args.mlst is None:
                msg.append('Please provide species name using --mlst')
    if args.minFrequencyDominantAllele < 0 or args.minFrequencyDominantAllele > 1:
        msg.append('--minFrequencyDominantAllele should be a value between [0, 1]')
    if args.minGeneCoverage < 0 or args.minGeneCoverage > 100:
        msg.append('--minGeneCoverage should be a value between [0, 100]')
    if args.minGeneIdentity < 0 or args.minGeneIdentity > 100:
        msg.append('--minGeneIdentity should be a value between [0, 100]')
    if args.notWriteConsensus and args.doubleRun:
        msg.append('--notWriteConsensus and --doubleRun cannot be used together.'
                   ' Maybe you only want to use --doubleRun')

    if len(msg) > 0:
        argparse.ArgumentParser.error('\n'.join(msg))

    start_time = time.time()

    number_samples_successfully, samples_total_number = run_rematch(args)

    print('\n' + 'END ReMatCh')
    print('\n' +
          str(number_samples_successfully) + ' samples out of ' + str(samples_total_number) + ' run successfully')
    time_taken = utils.run_time(start_time)
    del time_taken

    if number_samples_successfully == 0:
        sys.exit('No samples run successfully!')


if __name__ == "__main__":
    main()

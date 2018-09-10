import utils
import os.path
import multiprocessing
import sys
import functools
import time
import subprocess


def get_read_run_info(ena_id):
    import urllib

    url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=' + ena_id + '&result=read_run'

    read_run_info = None
    try:
        url = urllib.urlopen(url)
        read_run_info = url.read().splitlines()
        if len(read_run_info) <= 1:
            read_run_info = None
    except Exception as error:
        print(error)

    return read_run_info


def get_download_information(read_run_info):
    header_line = read_run_info[0].split('\t')
    info_line = read_run_info[1].split('\t')

    download_information = {'fastq': None, 'submitted': None, 'cram_index': None}
    download_types = ['aspera', 'ftp']

    for i in range(0, len(header_line)):
        header = header_line[i].lower().rsplit('_', 1)
        if header[0] in download_information.keys():
            if header[1] in download_types:
                if len(info_line[i]) > 0:
                    files_path = info_line[i].split(';')
                    if len(files_path) > 2:
                        print('WARNING: Were found more files than expected in'
                              ' {download_information}-{download_types} download'
                              ' links!'.format(download_information=header[0], download_types=header[1]))
                    if download_information[header[0]] is None:
                        download_information[header[0]] = {}
                    download_information[header[0]][header[1]] = files_path

    return download_information


def get_sequencing_information(read_run_info):
    header_line = read_run_info[0].split('\t')
    info_line = read_run_info[1].split('\t')

    sequencing_information = {'run_accession': None, 'instrument_platform': None, 'instrument_model': None,
                              'library_layout': None, 'library_source': None, 'extra_run_accession': None,
                              'nominal_length': None, 'read_count': None, 'base_count': None,
                              'date_download': time.strftime("%Y-%m-%d")}

    for i in range(0, len(header_line)):
        header = header_line[i].lower()
        if header in sequencing_information.keys():
            if len(info_line[i]) > 0:
                sequencing_information[header] = info_line[i]

    if len(read_run_info) > 2:
        extra_run_accession = []
        for i in range(2, len(read_run_info)):
            info = read_run_info[i].split('\t')
            for j in range(0, len(header_line)):
                header = header_line[j].lower()
                if header == 'run_accession':
                    if len(info[j]) > 0:
                        extra_run_accession.append(info[j])
        if len(extra_run_accession) >= 1:
            sequencing_information['extra_run_accession'] = ','.join(extra_run_accession)

    return sequencing_information


@utils.trace_unhandled_exceptions
def download_with_aspera(aspera_file_path, aspera_key, outdir, pickle_prefix, sra, ena_id):
    command = ['ascp', '-QT', '-l', '300m', '', '-i', aspera_key, '', outdir]
    if not sra:
        command[4] = '-P33001'
        command[7] = str('era-fasp@' + aspera_file_path)
        pickle = pickle_prefix + '.' + aspera_file_path.rsplit('/', 1)[1]
    else:
        command[7] = 'anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/{a}/{b}/{c}/{c}.sra'.format(
            a=ena_id[:3], b=ena_id[:6], c=ena_id)
        pickle = pickle_prefix + '.' + ena_id

    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, 3600, True)

    utils.saveVariableToPickle(run_successfully, outdir, pickle)


@utils.trace_unhandled_exceptions
def download_with_wget(ftp_file_path, outdir, pickle_prefix, sra, ena_id):
    command = ['wget', '--tries=1', '', '-O', '']
    if not sra:
        command[2] = ftp_file_path
        file_download = ftp_file_path.rsplit('/', 1)[1]
        command[4] = os.path.join(outdir, file_download)
        pickle = pickle_prefix + '.' + file_download
    else:
        command[2] = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{a}/{b}/{c}/{c}.sra'.format(
            a=ena_id[:3], b=ena_id[:6], c=ena_id)
        command[4] = os.path.join(outdir, ena_id + '.sra')
        pickle = pickle_prefix + '.' + ena_id
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, 3600, True)

    utils.saveVariableToPickle(run_successfully, outdir, pickle)


@utils.trace_unhandled_exceptions
def download_with_sra_prefetch(aspera_key, outdir, pickle_prefix, ena_id):
    command = ['prefetch', '', ena_id]

    if aspera_key is not None:
        _, ascp, _ = utils.runCommandPopenCommunicate(['which', 'ascp'], False, None, False)
        command[1] = '-a {ascp}|{aspera_key}'.format(ascp=ascp.splitlines()[0], aspera_key=aspera_key)

    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, 3600, True)
    if run_successfully:
        _, prefetch_outdir, _ = utils.runCommandPopenCommunicate(['echo', '$HOME/ncbi/public/sra'], True, None, False)

        try:
            os.rename(os.path.join(prefetch_outdir.splitlines()[0], ena_id + '.sra'),
                      os.path.join(outdir, ena_id + '.sra'))
        except OSError as e:
            print('Found the following error:'
                  '{}'.format(e))

            from shutil import copy as shutil_copy

            shutil_copy(os.path.join(prefetch_outdir.splitlines()[0], ena_id + '.sra'),
                        os.path.join(outdir, ena_id + '.sra'))
            os.remove(os.path.join(prefetch_outdir.splitlines()[0], ena_id + '.sra'))

    utils.saveVariableToPickle(run_successfully, outdir, pickle_prefix + '.' + ena_id)


@utils.trace_unhandled_exceptions
def download_with_curl(ftp_file_path, outdir, pickle_prefix, sra, ena_id):
    command = ['curl', '--retry', '1', '', '-o', '']
    if not sra:
        command[3] = ftp_file_path
        file_download = ftp_file_path.rsplit('/', 1)[1]
        command[5] = os.path.join(outdir, file_download)
        pickle = pickle_prefix + '.' + file_download
    else:
        command[3] = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{a}/{b}/{c}/{c}.sra'.format(
            a=ena_id[:3], b=ena_id[:6], c=ena_id)
        command[5] = os.path.join(outdir, ena_id + '.sra')
        pickle = pickle_prefix + '.' + ena_id
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, 3600, True)

    utils.saveVariableToPickle(run_successfully, outdir, pickle)


def get_pickle_run_successfully(directory, pickle_prefix):
    run_successfully = True
    read_pickle = False

    files = find_files(directory, pickle_prefix, '.pkl')
    if files is not None:
        for file_found in files:
            if run_successfully:
                run_successfully = utils.extractVariableFromPickle(file_found)
                read_pickle = True

            os.remove(file_found)

    if not read_pickle:
        run_successfully = False

    return run_successfully


def curl_installed():
    command = ['which', 'curl']
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
    return run_successfully


def download(download_information_type, aspera_key, outdir, sra, sra_opt, ena_id):
    pickle_prefix = 'download'

    run_successfully = False
    download_sra = False

    if not sra:
        if aspera_key is not None and download_information_type['aspera'] is not None:
            pool = multiprocessing.Pool(processes=2)
            for file_download in download_information_type['aspera']:
                pool.apply_async(download_with_aspera, args=(file_download, aspera_key, outdir, pickle_prefix, sra,
                                                             ena_id,))
            pool.close()
            pool.join()
            run_successfully = get_pickle_run_successfully(outdir, pickle_prefix)
        if not run_successfully and download_information_type['ftp'] is not None:
            if curl_installed():
                pool = multiprocessing.Pool(processes=2)
                for file_download in download_information_type['ftp']:
                    pool.apply_async(download_with_curl, args=(file_download, outdir, pickle_prefix, sra, ena_id,))
                pool.close()
                pool.join()
                run_successfully = get_pickle_run_successfully(outdir, pickle_prefix)
            if not run_successfully:
                pool = multiprocessing.Pool(processes=2)
                for file_download in download_information_type['ftp']:
                    pool.apply_async(download_with_wget, args=(file_download, outdir, pickle_prefix, sra, ena_id,))
                pool.close()
                pool.join()
                run_successfully = get_pickle_run_successfully(outdir, pickle_prefix)

    if not run_successfully and (sra or sra_opt):
        if aspera_key is not None:
            download_with_aspera(None, aspera_key, outdir, pickle_prefix, sra or sra_opt, ena_id)
            run_successfully = get_pickle_run_successfully(outdir, pickle_prefix)
        if not run_successfully:
            download_with_sra_prefetch(aspera_key, outdir, pickle_prefix, ena_id)
            run_successfully = get_pickle_run_successfully(outdir, pickle_prefix)
            if not run_successfully:
                if curl_installed():
                    download_with_curl(None, outdir, pickle_prefix, sra or sra_opt, ena_id)
                    run_successfully = get_pickle_run_successfully(outdir, pickle_prefix)
                if not run_successfully:
                    download_with_wget(None, outdir, pickle_prefix, sra or sra_opt, ena_id)
                    run_successfully = get_pickle_run_successfully(outdir, pickle_prefix)

        if run_successfully:
            download_sra = True

    return run_successfully, download_sra


def download_files(download_information, aspera_key, outdir, download_cram_bam_true, sra, sra_opt, ena_id):
    run_successfully = False
    cram_index_run_successfully = False
    download_sra = False

    if download_information['fastq'] is not None:
        run_successfully, download_sra = download(download_information['fastq'], aspera_key, outdir, sra, sra_opt,
                                                  ena_id)

    if not run_successfully:
        if download_information['submitted'] is not None:
            if not download_cram_bam_true:
                cram_bam = False
                for i in download_information['submitted']:
                    if download_information['submitted'][i][0].endswith(('.cram', '.bam')):
                        cram_bam = True
                        break
                if not cram_bam:
                    run_successfully, download_sra = download(download_information['submitted'], aspera_key, outdir,
                                                              False, False, ena_id)

            elif download_cram_bam_true:
                run_successfully, download_sra = download(download_information['submitted'], aspera_key, outdir, False,
                                                          False, ena_id)
                if run_successfully and download_information['cram_index'] is not None:
                    cram_index_run_successfully = download(download_information['cram_index'], aspera_key, outdir,
                                                           False, False, ena_id)

    if not run_successfully and (sra or sra_opt):
        run_successfully, download_sra = download(download_information['fastq'], aspera_key, outdir, True, sra_opt,
                                                  ena_id)

    return run_successfully, cram_index_run_successfully, download_sra


def sort_alignment(alignment_file, output_file, sort_by_name_true, threads):
    out_format_string = os.path.splitext(output_file)[1][1:].lower()
    command = ['samtools', 'sort', '-o', output_file, '-O', out_format_string, '', '-@', str(threads), alignment_file]
    if sort_by_name_true:
        command[6] = '-n'
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    if not run_successfully:
        output_file = None

    return run_successfully, output_file


def alignment_to_fastq(alignment_file, threads, pair_end_type):
    fastq_basename = os.path.splitext(alignment_file)[0]
    outfiles = None
    bam_file = fastq_basename + '.temp.bam'
    # sort cram
    run_successfully, bam_file = sort_alignment(alignment_file, bam_file, True, threads)
    if run_successfully:
        command = ['samtools', 'fastq', '', bam_file]
        if pair_end_type.lower() == 'paired':
            command[2] = '-1 ' + str(fastq_basename + '_1.fq') + ' -2 ' + str(fastq_basename + '_2.fq')
        elif pair_end_type == 'single':
            command[2] = '-0 ' + str(fastq_basename + '.fq')

        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
        if run_successfully:
            if pair_end_type.lower() == 'paired':
                outfiles = [str(fastq_basename + '_1.fq'), str(fastq_basename + '_2.fq')]
            elif pair_end_type.lower() == 'single':
                outfiles = [str(fastq_basename + '.fq')]

    if bam_file is not None and os.path.isfile(bam_file):
        os.remove(bam_file)

    return run_successfully, outfiles


def formart_fastq_headers(in_fastq_1, in_fastq_2):
    import itertools

    out_fastq_1 = in_fastq_1 + '.temp'
    out_fastq_2 = in_fastq_2 + '.temp'
    writer_in_fastq_1 = open(out_fastq_1, 'wt')
    writer_in_fastq_2 = open(out_fastq_2, 'wt')
    outfiles = [out_fastq_1, out_fastq_2]
    with open(in_fastq_1, 'rtU') as reader_in_fastq_1, open(in_fastq_2, 'rtU') as reader_in_fastq_2:
        plus_line = True
        quality_line = True
        number_reads = 0
        for in_1, in_2 in itertools.izip(reader_in_fastq_1, reader_in_fastq_2):
            if len(in_1) > 0:
                in_1 = in_1.splitlines()[0]
                in_2 = in_2.splitlines()[0]
                if in_1.startswith('@') and plus_line and quality_line:
                    if in_1 != in_2:
                        sys.exit('The PE fastq files are not aligned properly!')
                    in_1 += '/1' + '\n'
                    in_2 += '/2' + '\n'
                    writer_in_fastq_1.write(in_1)
                    writer_in_fastq_2.write(in_2)
                    plus_line = False
                    quality_line = False
                elif in_1.startswith('+') and not plus_line:
                    in_1 += '\n'
                    writer_in_fastq_1.write(in_1)
                    writer_in_fastq_2.write(in_1)
                    plus_line = True
                elif plus_line and not quality_line:
                    in_1 += '\n'
                    in_2 += '\n'
                    writer_in_fastq_1.write(in_1)
                    writer_in_fastq_2.write(in_2)
                    writer_in_fastq_1.flush()
                    writer_in_fastq_2.flush()
                    number_reads += 1
                    quality_line = True
                else:
                    in_1 += '\n'
                    in_2 += '\n'
                    writer_in_fastq_1.write(in_1)
                    writer_in_fastq_2.write(in_2)
    return number_reads, outfiles


@utils.trace_unhandled_exceptions
def gzip_files(file_2_compress, pickle_prefix, outdir):
    if file_2_compress.endswith('.temp'):
        out_file = os.path.splitext(file_2_compress)[0]
    else:
        out_file = file_2_compress

    command = ['gzip', '--stdout', '--best', file_2_compress, '>', str(out_file + '.gz')]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, True)
    if run_successfully:
        os.remove(file_2_compress)

    utils.saveVariableToPickle(run_successfully, outdir, str(pickle_prefix + '.' + os.path.basename(file_2_compress)))


def find_files(directory, prefix, suffix):
    list_files_found = []
    files = [f for f in os.listdir(directory) if not f.startswith('.') and os.path.isfile(os.path.join(directory, f))]
    for file_found in files:
        if file_found.startswith(prefix) and file_found.endswith(suffix):
            file_path = os.path.join(directory, file_found)
            list_files_found.append(file_path)

    if len(list_files_found) == 0:
        list_files_found = None

    return list_files_found


def compress_files(fastq_files, outdir, threads):
    pickle_prefix = 'compress'
    compressed_fastq_files = None

    pool = multiprocessing.Pool(processes=threads)
    for fastq in fastq_files:
        pool.apply_async(gzip_files, args=(fastq, pickle_prefix, outdir,))
    pool.close()
    pool.join()

    run_successfully = get_pickle_run_successfully(outdir, pickle_prefix)
    if run_successfully:
        compressed_fastq_files = find_files(outdir, '', '.gz')

    return run_successfully, compressed_fastq_files


def bam_cram_2_fastq(alignment_file, outdir, threads, pair_end_type):
    run_successfully, fastq_files = alignment_to_fastq(alignment_file, threads, pair_end_type)
    if run_successfully:
        if pair_end_type.lower() == 'paired':
            number_reads, fastq_files = formart_fastq_headers(fastq_files[0], fastq_files[1])

        run_successfully, fastq_files = compress_files(fastq_files, outdir, threads)

    return run_successfully, fastq_files


def check_correct_links(download_information):
    for i in download_information:
        if download_information[i] is not None:
            if download_information[i]['aspera'] is not None:
                for j in range(0, len(download_information[i]['aspera'])):
                    if download_information[i]['aspera'][j].startswith('fasp.sra.ebi.ac.uk/'):
                        download_information[i]['aspera'][j] = download_information[i]['aspera'][j].replace(
                            'fasp.sra.ebi.ac.uk/', 'fasp.sra.ebi.ac.uk:/', 1)
            if download_information[i]['ftp'] is not None:
                for j in range(0, len(download_information[i]['ftp'])):
                    if '#' in download_information[i]['ftp'][j]:
                        download_information[i]['ftp'][j] = download_information[i]['ftp'][j].replace('#', '%23')
    return download_information


def get_fastq_files(download_dir, cram_index_run_successfully, threads, download_paired_type):
    run_successfully = False
    downloaded_files = find_files(download_dir, '', '')
    if cram_index_run_successfully:
        cram_file = None
        for i in downloaded_files:
            if i.endswith('.cram'):
                cram_file = i
        run_successfully, downloaded_files = bam_cram_2_fastq(cram_file, download_dir, threads, download_paired_type)
    else:
        if downloaded_files is not None and len(downloaded_files) > 0:
            run_successfully = True

    return run_successfully, downloaded_files


def rename_move_files(list_files, new_name, outdir, download_paired_type):
    list_new_files = {}
    run_successfully = False

    for i in range(0, len(list_files)):
        temp_name = utils.rchop(os.path.basename(list_files[i]), 'astq.gz')
        if len(temp_name) == len(os.path.basename(list_files[i])):
            temp_name = utils.rchop(os.path.basename(list_files[i]), 'q.gz')
        if download_paired_type.lower() == 'paired':
            if temp_name.endswith(('_R1_001.f', '_1.f')):
                list_new_files[i] = os.path.join(outdir, new_name + '_1.fq.gz')
            elif temp_name.endswith(('_R2_001.f', '_2.f')):
                list_new_files[i] = os.path.join(outdir, new_name + '_2.fq.gz')
        else:
            if not temp_name.endswith(('_R1_001.f', '_R2_001.f')):
                list_new_files[i] = os.path.join(outdir, new_name + '.fq.gz')
                if temp_name.endswith(('_1.f', '_2.f')):
                    print('WARNING: possible single-end file conflict with pair-end (' + list_files[i] + ')!')

    if len(list_new_files) == 2 and download_paired_type.lower() == 'paired':
        run_successfully = True
    elif len(list_new_files) == 1 and download_paired_type.lower() == 'single':
        run_successfully = True

    if run_successfully:
        try:
            for i in range(0, len(list_files)):
                if i not in list_new_files:
                    if os.path.isfile(list_files[i]):
                        os.remove(list_files[i])
                else:
                    os.rename(list_files[i], list_new_files[i])
            list_new_files = list_new_files.values()
        except Exception as e:
            print(e)
            run_successfully = False

    if not run_successfully:
        list_new_files = None

    return run_successfully, list_new_files


# @utils.trace_unhandled_exceptions
def rename_header_sra(fastq):
    run_successfully = False
    try:
        command = ['gawk', '\'{if(NR%4==1) $0=gensub(/\./, \"/\", 2); print}\'', fastq, '|', 'gzip', '-1', '>',
                   str(fastq + '.gz')]
        print('Running: ' + str(' '.join(command)))
        return_code = subprocess.call(' '.join(command), shell=True)
        if return_code == 0:
            run_successfully = True
        else:
            print('Something went wrong with command: {command}'.format(command=' '.join(command)))
    except Exception as e:
        print(e)

    return run_successfully


def sra_2_fastq(download_dir, ena_id):
    command = ['fastq-dump', '-I', '-O', download_dir, '--split-files', '{download_dir}{ena_id}.sra'.format(
        download_dir=download_dir, ena_id=ena_id)]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, 3600, True)
    if run_successfully:
        files = [os.path.join(download_dir, f) for f in os.listdir(download_dir)
                 if not f.startswith('.') and os.path.isfile(os.path.join(download_dir, f)) and f.endswith('.fastq')]

        pool = multiprocessing.Pool(processes=2)
        results = []
        p = pool.map_async(rename_header_sra, files, callback=results.extend)
        p.wait()

        run_successfully = all(results)

    return run_successfully


download_timer = functools.partial(utils.timer, name='Download module')


@download_timer
def run_download(ena_id, download_paired_type, aspera_key, outdir, download_cram_bam_true, threads, instrument_platform,
                 sra, sra_opt):
    download_dir = os.path.join(outdir, 'download', '')
    utils.removeDirectory(download_dir)
    os.mkdir(download_dir)

    run_successfully = False
    downloaded_files = None
    sequencing_information = {'run_accession': None, 'instrument_platform': None, 'instrument_model': None,
                              'library_layout': None, 'library_source': None, 'extra_run_accession': None,
                              'nominal_length': None, 'read_count': None, 'base_count': None,
                              'date_download': time.strftime("%Y-%m-%d")}

    read_run_info = get_read_run_info(ena_id)
    if read_run_info is not None:
        download_information = get_download_information(read_run_info)
        download_information = check_correct_links(download_information)
        sequencing_information = get_sequencing_information(read_run_info)

        if instrument_platform.lower() == 'all' or \
                (sequencing_information['instrument_platform'] is not None and
                 sequencing_information['instrument_platform'].lower() == instrument_platform.lower()):
            if download_paired_type.lower() == 'both' or \
                    (sequencing_information['library_layout'] is not None and
                     sequencing_information['library_layout'].lower() == download_paired_type.lower()):
                run_successfully, cram_index_run_successfully, download_sra = download_files(download_information,
                                                                                             aspera_key, download_dir,
                                                                                             download_cram_bam_true,
                                                                                             sra, sra_opt, ena_id)
                if download_sra:
                    run_successfully = sra_2_fastq(download_dir, ena_id)
                if run_successfully:
                    run_successfully, downloaded_files = get_fastq_files(download_dir, cram_index_run_successfully,
                                                                         threads,
                                                                         sequencing_information['library_layout'])
                if run_successfully and downloaded_files is not None:
                    run_successfully, downloaded_files = rename_move_files(downloaded_files,
                                                                           sequencing_information['run_accession'],
                                                                           outdir,
                                                                           sequencing_information['library_layout'])
    else:
        if sra or sra_opt:
            run_successfully, cram_index_run_successfully, download_sra = download_files({'fastq': None,
                                                                                          'submitted': None,
                                                                                          'cram_index': None},
                                                                                         aspera_key, download_dir,
                                                                                         download_cram_bam_true, sra,
                                                                                         sra_opt, ena_id)
            if download_sra:
                run_successfully = sra_2_fastq(download_dir, ena_id)
            if run_successfully:
                run_successfully, downloaded_files = get_fastq_files(download_dir, cram_index_run_successfully, threads,
                                                                     'paired')
                if not run_successfully:
                    run_successfully, downloaded_files = get_fastq_files(download_dir, cram_index_run_successfully,
                                                                         threads, 'single')
            if run_successfully and downloaded_files is not None:
                run_successfully, downloaded_files = rename_move_files(downloaded_files, ena_id, outdir, 'paired')
                if not run_successfully:
                    run_successfully, downloaded_files = rename_move_files(downloaded_files, ena_id, outdir, 'single')

    utils.removeDirectory(download_dir)

    return run_successfully, downloaded_files, sequencing_information

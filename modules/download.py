import utils
import os.path
import multiprocessing
import sys
import functools
import time


def getReadRunInfo(ena_id):
    import urllib

    url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=' + ena_id + '&result=read_run'

    readRunInfo = None
    try:
        url = urllib.urlopen(url)
        readRunInfo = url.read().splitlines()
        if len(readRunInfo) <= 1:
            readRunInfo = None
    except Exception as error:
        print error

    return readRunInfo


def getDownloadInformation(readRunInfo):
    header_line = readRunInfo[0].split('\t')
    info_line = readRunInfo[1].split('\t')

    downloadInformation = {'fastq': None, 'submitted': None, 'cram_index': None}
    download_types = ['aspera', 'ftp']

    for i in range(0, len(header_line)):
        header = header_line[i].lower().rsplit('_', 1)
        if header[0] in downloadInformation.keys():
            if header[1] in download_types:
                if len(info_line[i]) > 0:
                    files_path = info_line[i].split(';')
                    if len(files_path) > 2:
                        print 'WARNING: Were found more files than expected in ' + header[1] + ' download links!'
                    if downloadInformation[header[0]] is None:
                        downloadInformation[header[0]] = {}
                    downloadInformation[header[0]][header[1]] = files_path

    return downloadInformation


def getSequencingInformation(readRunInfo):
    header_line = readRunInfo[0].split('\t')
    info_line = readRunInfo[1].split('\t')

    sequencingInformation = {'run_accession': None, 'instrument_platform': None, 'instrument_model': None, 'library_layout': None, 'library_source': None, 'extra_run_accession': None, 'nominal_length': None, 'read_count': None, 'base_count': None}

    for i in range(0, len(header_line)):
        header = header_line[i].lower()
        if header in sequencingInformation.keys():
            if len(info_line[i]) > 0:
                sequencingInformation[header] = info_line[i]

    if len(readRunInfo) > 2:
        extra_run_accession = []
        for i in range(2, len(readRunInfo)):
            info = readRunInfo[i].split('\t')
            for j in range(0, len(header_line)):
                header = header_line[j].lower()
                if header == 'run_accession':
                    if len(info[j]) > 0:
                        extra_run_accession.append(info[j])
        if len(extra_run_accession) >= 1:
            sequencingInformation['extra_run_accession'] = ','.join(extra_run_accession)

    return sequencingInformation


@utils.trace_unhandled_exceptions
def downloadWithAspera(aspera_file_path, asperaKey, outdir, pickle_prefix):
    command = ['ascp', '-QT', '-l', '300m', '-P33001', '-i', asperaKey, str('era-fasp@' + aspera_file_path), outdir]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, 3600, True)

    utils.saveVariableToPickle(run_successfully, outdir, str(pickle_prefix + '.' + aspera_file_path.rsplit('/', 1)[1]))


@utils.trace_unhandled_exceptions
def downloadWithWget(ftp_file_path, outdir, pickle_prefix):
    file_download = ftp_file_path.rsplit('/', 1)[1]
    command = ['wget', '--tries=2', ftp_file_path, '-O', os.path.join(outdir, file_download)]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, 3600, True)

    utils.saveVariableToPickle(run_successfully, outdir, str(pickle_prefix + '.' + file_download))


@utils.trace_unhandled_exceptions
def downloadWithCurl(ftp_file_path, outdir, pickle_prefix):
    file_download = ftp_file_path.rsplit('/', 1)[1]
    command = ['curl', '--retry', '2', ftp_file_path, '-O', os.path.join(outdir, file_download)]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, 3600, True)

    utils.saveVariableToPickle(run_successfully, outdir, str(pickle_prefix + '.' + file_download))


def getPickleRunSuccessfully(directory, pickle_prefix):
    run_successfully = True
    read_pickle = False

    files = findFiles(directory, pickle_prefix, '.pkl')
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


def download(downloadInformation_type, asperaKey, outdir):
    pickle_prefix = 'download'

    run_successfully = False

    if asperaKey is not None and downloadInformation_type['aspera'] is not None:
        pool = multiprocessing.Pool(processes=2)
        for file_download in downloadInformation_type['aspera']:
            pool.apply_async(downloadWithAspera, args=(file_download, asperaKey, outdir, pickle_prefix,))
        pool.close()
        pool.join()

        run_successfully = getPickleRunSuccessfully(outdir, pickle_prefix)

    if downloadInformation_type['ftp'] is not None and not run_successfully:
        if curl_installed():
            pool = multiprocessing.Pool(processes=2)
            for file_download in downloadInformation_type['ftp']:
                pool.apply_async(downloadWithCurl, args=(file_download, outdir, pickle_prefix,))
            pool.close()
            pool.join()
            run_successfully = getPickleRunSuccessfully(outdir, pickle_prefix)
        if not run_successfully:
            pool = multiprocessing.Pool(processes=2)
            for file_download in downloadInformation_type['ftp']:
                pool.apply_async(downloadWithWget, args=(file_download, outdir, pickle_prefix,))
            pool.close()
            pool.join()
            run_successfully = getPickleRunSuccessfully(outdir, pickle_prefix)

    return run_successfully


def downloadFiles(downloadInformation, asperaKey, outdir, download_cram_bam_True):
    run_successfully = False
    cram_index_run_successfully = False

    if downloadInformation['fastq'] is not None:
        run_successfully = download(downloadInformation['fastq'], asperaKey, outdir)

    if not run_successfully:
        if downloadInformation['submitted'] is not None:
            if not download_cram_bam_True:
                cram_bam = False
                for i in downloadInformation['submitted']:
                    if downloadInformation['submitted'][i][0].endswith(('.cram', '.bam')):
                        cram_bam = True
                        break
                if not cram_bam:
                    run_successfully = download(downloadInformation['submitted'], asperaKey, outdir)

            elif download_cram_bam_True:
                run_successfully = download(downloadInformation['submitted'], asperaKey, outdir)
                if run_successfully and downloadInformation['cram_index'] is not None:
                    cram_index_run_successfully = download(downloadInformation['cram_index'], asperaKey, outdir)

    return run_successfully, cram_index_run_successfully


def sortAlignment(alignment_file, output_file, sortByName_True, threads):
    outFormat_string = os.path.splitext(output_file)[1][1:].lower()
    command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
    if sortByName_True:
        command[6] = '-n'
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    if not run_successfully:
        output_file = None

    return run_successfully, output_file


def alignmentToFastq(alignment_file, outdir, threads, pair_end_type):
    fastq_basename = os.path.splitext(alignment_file)[0]
    outfiles = None
    bamFile = fastq_basename + '.temp.bam'
    # sort cram
    run_successfully, bamFile = sortAlignment(alignment_file, bamFile, True, threads)
    if run_successfully:
        command = ['samtools', 'fastq', '', bamFile]
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

    if os.path.isfile(bamFile):
        os.remove(bamFile)

    return run_successfully, outfiles


def formartFastqHeaders(in_fastq_1, in_fastq_2):
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
def gzipFiles(file_2_compress, pickle_prefix, outdir):
    out_file = None
    if file_2_compress.endswith('.temp'):
        out_file = os.path.splitext(file_2_compress)[0]
    else:
        out_file = file_2_compress

    command = ['gzip', '--stdout', '--best', file_2_compress, '>', str(out_file + '.gz')]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, True)
    if run_successfully:
        os.remove(file_2_compress)

    utils.saveVariableToPickle(run_successfully, outdir, str(pickle_prefix + '.' + os.path.basename(file_2_compress)))


def findFiles(directory, prefix, suffix):
    list_files_found = []
    files = [f for f in os.listdir(directory) if not f.startswith('.') and os.path.isfile(os.path.join(directory, f))]
    for file_found in files:
        if file_found.startswith(prefix) and file_found.endswith(suffix):
            file_path = os.path.join(directory, file_found)
            list_files_found.append(file_path)

    if len(list_files_found) == 0:
        list_files_found = None

    return list_files_found


def compressFiles(fastq_files, outdir, threads):
    pickle_prefix = 'compress'
    compressed_fastq_files = None

    pool = multiprocessing.Pool(processes=threads)
    for fastq in fastq_files:
        pool.apply_async(gzipFiles, args=(fastq, pickle_prefix, outdir,))
    pool.close()
    pool.join()

    run_successfully = getPickleRunSuccessfully(outdir, pickle_prefix)
    if run_successfully:
        compressed_fastq_files = findFiles(outdir, '', '.gz')

    return run_successfully, compressed_fastq_files


def bamCram_2_fastq(alignment_file, outdir, threads, pair_end_type):
    run_successfully, fastq_files = alignmentToFastq(alignment_file, outdir, threads, pair_end_type)
    if run_successfully:
        if pair_end_type.lower() == 'paired':
            number_reads, fastq_files = formartFastqHeaders(fastq_files[0], fastq_files[1])

        run_successfully, fastq_files = compressFiles(fastq_files, outdir, threads)

    return run_successfully, fastq_files


def check_correct_links(downloadInformation):
    for i in downloadInformation:
        if downloadInformation[i] is not None:
            if downloadInformation[i]['aspera'] is not None:
                for j in range(0, len(downloadInformation[i]['aspera'])):
                    if downloadInformation[i]['aspera'][j].startswith('fasp.sra.ebi.ac.uk/'):
                        downloadInformation[i]['aspera'][j] = downloadInformation[i]['aspera'][j].replace('fasp.sra.ebi.ac.uk/', 'fasp.sra.ebi.ac.uk:/', 1)
            if downloadInformation[i]['ftp'] is not None:
                for j in range(0, len(downloadInformation[i]['ftp'])):
                    if '#' in downloadInformation[i]['ftp'][j]:
                        downloadInformation[i]['ftp'][j] = downloadInformation[i]['ftp'][j].replace('#', '%23')
    return downloadInformation


def get_fastq_files(download_dir, cram_index_run_successfully, threads, download_paired_type):
    run_successfully = False
    downloaded_files = findFiles(download_dir, '', '')
    if cram_index_run_successfully:
        cram_file = None
        for i in downloaded_files:
            if i.endswith('.cram'):
                cram_file = i
        run_successfully, downloaded_files = bamCram_2_fastq(cram_file, download_dir, threads, download_paired_type)
    else:
        if len(downloaded_files) > 0:
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
                    print 'WARNING: possible single-end file conflict with pair-end (' + list_files[i] + ')!'

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
            print e
            run_successfully = False

    if not run_successfully:
        list_new_files = None

    return run_successfully, list_new_files


download_timer = functools.partial(utils.timer, name='Download module')


@download_timer
def runDownload(ena_id, download_paired_type, asperaKey, outdir, download_cram_bam_True, threads, instrument_platform):
    download_dir = os.path.join(outdir, 'download', '')
    utils.removeDirectory(download_dir)
    os.mkdir(download_dir)

    run_successfully = False
    downloaded_files = None
    sequencingInformation = {'run_accession': None, 'instrument_platform': None, 'instrument_model': None, 'library_layout': None, 'library_source': None, 'extra_run_accession': None, 'nominal_length': None, 'read_count': None, 'base_count': None, 'date_download': None}

    readRunInfo = getReadRunInfo(ena_id)
    if readRunInfo is not None:
        downloadInformation = getDownloadInformation(readRunInfo)
        downloadInformation = check_correct_links(downloadInformation)
        sequencingInformation = getSequencingInformation(readRunInfo)
        sequencingInformation['date_download'] = time.strftime("%Y-%m-%d")

        if instrument_platform.lower() == 'all' or (sequencingInformation['instrument_platform'] is not None and sequencingInformation['instrument_platform'].lower() == instrument_platform.lower()):
            if download_paired_type.lower() == 'both' or (sequencingInformation['library_layout'] is not None and sequencingInformation['library_layout'].lower() == download_paired_type.lower()):
                run_successfully, cram_index_run_successfully = downloadFiles(downloadInformation, asperaKey, download_dir, download_cram_bam_True)
                if run_successfully:
                    run_successfully, downloaded_files = get_fastq_files(download_dir, cram_index_run_successfully, threads, sequencingInformation['library_layout'])
                if run_successfully and downloaded_files is not None:
                    run_successfully, downloaded_files = rename_move_files(downloaded_files, sequencingInformation['run_accession'], outdir, sequencingInformation['library_layout'])

    utils.removeDirectory(download_dir)

    return run_successfully, downloaded_files, sequencingInformation

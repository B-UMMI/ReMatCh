import pickle
from traceback import format_exception as traceback_format_exception
import shlex
import subprocess
from threading import Timer
import shutil
import time
from functools import wraps as functools_wraps
import os.path
import sys


def start_logger(workdir):
    time_str = time.strftime("%Y%m%d-%H%M%S")
    sys.stdout = Logger(workdir, time_str)
    logfile = sys.stdout.getLogFile()
    return logfile, time_str


class Logger(object):
    def __init__(self, out_directory, time_str):
        self.logfile = os.path.join(out_directory, str('run.' + time_str + '.log'))
        self.terminal = sys.stdout
        self.log = open(self.logfile, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        pass

    def getLogFile(self):
        return self.logfile


def get_cpu_information(outdir, time_str):
    with open(os.path.join(outdir, 'cpu_information.' + time_str + '.cpu.txt'), 'wt') as writer:
        command = ['cat', '/proc/cpuinfo']
        run_successfully, stdout, stderr = run_command_popen_communicate(command, False, None, False)
        if run_successfully:
            writer.write(stdout)

    with open(os.path.join(outdir, 'cpu_information.' + time_str + '.slurm.txt'), 'wt') as writer:
        for environment in sorted(os.environ):
            if environment.startswith('SLURM_'):
                writer.write('#' + environment + '\n' + os.environ[environment] + '\n')


def setPATHvariable(doNotUseProvidedSoftware, script_path):
    path_variable = os.environ['PATH']
    script_folder = os.path.dirname(script_path)
    # Set path to use provided softwares
    if not doNotUseProvidedSoftware:
        bowtie2 = os.path.join(script_folder, 'src', 'bowtie2-2.2.9')
        samtools = os.path.join(script_folder, 'src', 'samtools-1.3.1', 'bin')
        bcftools = os.path.join(script_folder, 'src', 'bcftools-1.3.1', 'bin')

        os.environ['PATH'] = str(':'.join([bowtie2, samtools, bcftools, path_variable]))

    # Print PATH variable
    print('\n' + 'PATH variable:')
    print(os.environ['PATH'])


def checkPrograms(programs_version_dictionary):
    print('\n' + 'Checking dependencies...')
    programs = programs_version_dictionary
    which_program = ['which', '']
    listMissings = []
    for program in programs:
        which_program[1] = program
        run_successfully, stdout, stderr = run_command_popen_communicate(which_program, False, None, False)
        if not run_successfully:
            listMissings.append(program + ' not found in PATH.')
        else:
            print(stdout.splitlines()[0])
            if programs[program][0] is None:
                print(program + ' (impossible to determine programme version) found at: ' + stdout.splitlines()[0])
            else:
                if program.endswith('.jar'):
                    check_version = ['java', '-jar', stdout.splitlines()[0], programs[program][0]]
                    programs[program].append(stdout.splitlines()[0])
                else:
                    check_version = [stdout.splitlines()[0], programs[program][0]]
                run_successfully, stdout, stderr = run_command_popen_communicate(check_version, False, None, False)
                if stdout == '':
                    stdout = stderr
                if program in ['wget', 'awk']:
                    version_line = stdout.splitlines()[0].split(' ', 3)[2]
                elif program in ['prefetch', 'fastq-dump']:
                    version_line = stdout.splitlines()[1].split(' ')[-1]
                else:
                    version_line = stdout.splitlines()[0].split(' ')[-1]
                replace_characters = ['"', 'v', 'V', '+', ',']
                for i in replace_characters:
                    version_line = version_line.replace(i, '')
                print(program + ' (' + version_line + ') found')
                if programs[program][1] == '>=':
                    program_found_version = version_line.split('.')
                    program_version_required = programs[program][2].split('.')
                    if len(program_version_required) == 3:
                        if len(program_found_version) == 2:
                            program_found_version.append(0)
                        else:
                            program_found_version[2] = program_found_version[2].split('_')[0]
                    for i in range(0, len(program_version_required)):
                        if int(program_found_version[i]) > int(program_version_required[i]):
                            break
                        elif int(program_found_version[i]) == int(program_version_required[i]):
                            continue
                        else:
                            listMissings.append('It is required ' + program + ' with version ' +
                                                programs[program][1] + ' ' + programs[program][2])
                else:
                    if version_line != programs[program][2]:
                        listMissings.append('It is required ' + program + ' with version ' + programs[program][1] +
                                            ' ' + programs[program][2])
    return listMissings


def requiredPrograms(asperaKey, downloadCramBam, SRA, SRAopt):
    programs_version_dictionary = {}
    programs_version_dictionary['wget'] = ['--version', '>=', '1.12']
    programs_version_dictionary['gzip'] = ['--version', '>=', '1.6']
    programs_version_dictionary['bowtie2'] = ['--version', '>=', '2.2.9']
    programs_version_dictionary['samtools'] = ['--version', '==', '1.3.1']
    programs_version_dictionary['bcftools'] = ['--version', '==', '1.3.1']
    if asperaKey is not None:
        programs_version_dictionary['ascp'] = ['--version', '>=', '3.6.1']
    if SRA or SRAopt:
        programs_version_dictionary['prefetch'] = ['--version', '>=', '2.8.2']
        programs_version_dictionary['fastq-dump'] = ['--version', '>=', '2.8.2']
        programs_version_dictionary['awk'] = ['--version', '>=', '3.0.4']
    missingPrograms = checkPrograms(programs_version_dictionary)
    if len(missingPrograms) > 0:
        sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missingPrograms))


def general_information(logfile, version, outdir, time_str, doNotUseProvidedSoftware, asperaKey, downloadCramBam, SRA, SRAopt):
    # Check if output directory exists

    print('\n' + '==========> ReMatCh <==========')
    print('\n' + 'Program start: ' + time.ctime())

    # Tells where the logfile will be stored
    print('\n' + 'LOGFILE:')
    print(logfile)

    # Print command
    print('\n' + 'COMMAND:')
    script_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'rematch.py')
    print(sys.executable + ' ' + ' '.join(sys.argv))

    # Print directory where programme was lunch
    print('\n' + 'PRESENT DIRECTORY:')
    present_directory = os.path.abspath(os.getcwd())
    print(present_directory)

    # Print program version
    print('\n' + 'VERSION:')
    script_version_git(version, present_directory, script_path)

    # Get CPU information
    get_cpu_information(outdir, time_str)

    # Set and print PATH variable
    setPATHvariable(doNotUseProvidedSoftware, script_path)

    # Check programms
    requiredPrograms(asperaKey, downloadCramBam, SRA, SRAopt)

    return script_path


def script_version_git(version, current_directory, script_path, no_git_info=False):
    """
    Print script version and get GitHub commit information

    Parameters
    ----------
    version : str
        Version of the script, e.g. "4.0"
    current_directory : str
        Path to the directory where the script was start to run
    script_path : str
        Path to the script running
    no_git_info : bool, default False
        True if it is not necessary to retreive the GitHub commit information

    Returns
    -------

    """
    print('Version {}'.format(version))

    if not no_git_info:
        try:
            os.chdir(os.path.dirname(os.path.dirname(script_path)))
            command = ['git', 'log', '-1', '--date=local', '--pretty=format:"%h (%H) - Commit by %cn, %cd) : %s"']
            run_successfully, stdout, stderr = run_command_popen_communicate(command, False, 15, False)
            print(stdout)
            command = ['git', 'remote', 'show', 'origin']
            run_successfully, stdout, stderr = run_command_popen_communicate(command, False, 15, False)
            print(stdout)
        except:
            print('HARMLESS WARNING: git command possibly not found. The GitHub repository information will not be'
                  ' obtained.')
        finally:
            os.chdir(current_directory)


def run_time(start_time):
    end_time = time.time()
    time_taken = end_time - start_time
    hours, rest = divmod(time_taken, 3600)
    minutes, seconds = divmod(rest, 60)
    print('Runtime :' + str(hours) + 'h:' + str(minutes) + 'm:' + str(round(seconds, 2)) + 's')
    return round(time_taken, 2)


def timer(function, name):
    @functools_wraps(function)
    def wrapper(*args, **kwargs):
        print('\n' + 'RUNNING {0}\n'.format(name))
        start_time = time.time()

        results = list(function(*args, **kwargs))  # guarantees return is a list to allow .insert()

        time_taken = run_time(start_time)
        print('END {0}'.format(name))

        results.insert(0, time_taken)
        return results
    return wrapper


def remove_directory(directory):
    if os.path.isdir(directory):
        shutil.rmtree(directory)


def save_variable_to_pickle(variableToStore, outdir, prefix):
    pickleFile = os.path.join(outdir, str(prefix + '.pkl'))
    with open(pickleFile, 'wb') as writer:
        pickle.dump(variableToStore, writer)


def extract_variable_from_pickle(pickleFile):
    with open(pickleFile, 'rb') as reader:
        variable = pickle.load(reader)
    return variable


def trace_unhandled_exceptions(func):
    @functools_wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            print('Exception in ' + func.__name__)
            print(e)

            exc_type, exc_value, exc_tb = sys.exc_info()
            print(''.join(traceback_format_exception(exc_type, exc_value, exc_tb)))

            raise exc_type(exc_value)

    return wrapped_func


def kill_subprocess_Popen(subprocess_Popen, command):
    print('Command run out of time: ' + str(command))
    subprocess_Popen.kill()


def run_command_popen_communicate(command, shell_True, timeout_sec_None, print_comand_True):
    run_successfully = False
    if not isinstance(command, str):
        command = ' '.join(command)
    command = shlex.split(command)

    if print_comand_True:
        print('Running: ' + ' '.join(command))

    if shell_True:
        command = ' '.join(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    else:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    not_killed_by_timer = True
    if timeout_sec_None is None:
        stdout, stderr = proc.communicate()
    else:
        time_counter = Timer(timeout_sec_None, kill_subprocess_Popen, args=(proc, command,))
        time_counter.start()
        stdout, stderr = proc.communicate()
        time_counter.cancel()
        not_killed_by_timer = time_counter.isAlive()

    if proc.returncode == 0:
        run_successfully = True
    else:
        if not print_comand_True and not_killed_by_timer:
            print('Running: ' + str(command))
        if len(stdout) > 0:
            print('STDOUT')
            print(stdout.decode("utf-8"))
        if len(stderr) > 0:
            print('STDERR')
            print(stderr.decode("utf-8"))
    return run_successfully, stdout.decode("utf-8"), stderr.decode("utf-8")


def rchop(string, ending):
    if string.endswith(ending):
        string = string[:-len(ending)]
    return string


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    reverse_complement_string = ''

    seq = reversed(list(seq.upper()))

    for base in seq:
        reverse_complement_string += complement[base]

    return reverse_complement_string

import sys
import os
import urllib.request
import csv
from glob import glob
import re
import functools
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from . import utils
from . import rematch_module


def determine_species(species):
    species = species.lower().split(' ')

    if len(species) >= 2:
        species = species[:2]
        if species[1] in ('spp', 'spp.', 'complex'):
            species = [species[0]]

    return species


def check_existing_schema(species, schema_number, script_path):
    species = determine_species(species)

    if schema_number is None:
        schema_number = ''
    else:
        schema_number = '#' + str(schema_number)

    mlst_schemas_folder = os.path.join(os.path.dirname(script_path), 'modules', 'mlst_schemas', '')
    reference = []
    files = [f for f in os.listdir(mlst_schemas_folder) if not f.startswith('.') and os.path.isfile(os.path.join(mlst_schemas_folder, f))]
    for file_found in files:
        file_path = os.path.join(mlst_schemas_folder, file_found)
        if file_found.startswith('_'.join(species) + schema_number) and file_found.endswith('.fasta'):
            reference = file_path

    if len(reference) > 1:
        if schema_number == '':
            schema_number = '#1'
        for scheme in reference:
            if os.path.splitext(scheme)[0].endswith(schema_number):
                reference = [scheme]
                break
    if len(reference) == 0:
        reference = None
    elif len(reference) == 1:
        reference = reference[0]
    return reference


def write_mlst_reference(species, mlst_sequences, outdir, time_str):
    print('Writing MLST alleles as reference_sequences' + '\n')
    reference_file = os.path.join(outdir, str(species.replace('/', '_').replace(' ', '_') + '.' + time_str + '.fasta'))
    with open(reference_file, 'wt') as writer:
        for header, sequence in list(mlst_sequences.items()):
            writer.write('>' + header + '\n')
            fasta_sequence_lines = rematch_module.chunkstring(sequence, 80)
            for line in fasta_sequence_lines:
                writer.write(line + '\n')
    return reference_file


def get_st(mlst_dicts, dict_sequences):
    SequenceDict = mlst_dicts[0]
    STdict = mlst_dicts[1]
    lociOrder = mlst_dicts[2]

    alleles_profile = ['-'] * len(lociOrder)
    for x, sequence_data in list(dict_sequences.items()):
        if sequence_data['header'] not in SequenceDict:
            print(sequence_data['header'] + ' not found between consensus sequences!')
            break
        if sequence_data['sequence'] in list(SequenceDict[sequence_data['header']].keys()):
            allele_number = SequenceDict[sequence_data['header']][sequence_data['sequence']]
            alleles_profile[lociOrder.index(sequence_data['header'])] = allele_number
        else:
            for sequence_st, allele_number in list(SequenceDict[sequence_data['header']].items()):
                if sequence_st in sequence_data['sequence']:
                    alleles_profile[lociOrder.index(sequence_data['header'])] = allele_number

    alleles_profile = ','.join(alleles_profile)
    st = '-'
    if alleles_profile in STdict:
        st = STdict[alleles_profile]

    return st, alleles_profile


downloadPubMLST = functools.partial(utils.timer, name='Download PubMLST module')


@downloadPubMLST
def download_pub_mlst_xml(originalSpecies, schema_number, outdir):
    print('Searching MLST database for ' + originalSpecies)

    xmlURL = 'http://pubmlst.org/data/dbases.xml'
    try:
        content = urllib.request.urlopen(xmlURL)
        xml = content.read()
        tree = ET.fromstring(xml)
    except:
        print("Ooops! There might be a problem with the PubMLST service, try later or check if the xml is well formated"
              " at " + xmlURL)
        raise

    xmlData = {}

    if schema_number is None:
        schema_number = 1

    success = 0
    for scheme in tree.findall('species'):
        species_scheme = scheme.text.rstrip('\r\n').rsplit('#', 1)
        number_scheme = species_scheme[1] if len(species_scheme) == 2 else 1
        species_scheme = species_scheme[0]
        if determine_species(species_scheme) == determine_species(originalSpecies):
            if schema_number == number_scheme:
                success += 1
                xmlData[scheme.text.strip()] = {}
                for info in scheme:  # mlst
                    for database in info:  # database
                        for retrievedDate in database.findall('retrieved'):
                            retrieved = retrievedDate.text
                            xmlData[scheme.text.strip()][retrieved] = []
                        for profile in database.findall('profiles'):
                                profileURl = profile.find('url').text
                                xmlData[scheme.text.strip()][retrieved].append(profileURl)
                        for lociScheme in database.findall('loci'):
                            loci = {}
                            for locus in lociScheme:
                                locusID = locus.text
                                for locusInfo in locus:
                                    locusUrl = locusInfo.text
                                    loci[locusID.strip()] = locusUrl
                                xmlData[scheme.text.strip()][retrieved].append(loci)
    if success == 0:
        sys.exit("\tError. No schema found for %s. Please refer to https://pubmlst.org/databases/" % (originalSpecies))
    elif success > 1:
        keys = list(xmlData.keys())
        keys = sorted(keys)
        print("\tWarning. More than one schema found for %s. only keeping the first"
              " one... %s" % (originalSpecies, keys[0]))
        for key in keys[1:]:
            del xmlData[key]

    pubmlst_dir = os.path.join(outdir, 'pubmlst', '')
    if not os.path.isdir(pubmlst_dir):
        os.makedirs(pubmlst_dir)

    for SchemaName, info in list(xmlData.items()):
        STdict = {}
        SequenceDict = {}
        mlst_sequences = {}

        species_name = '_'.join(determine_species(SchemaName)).replace('/', '_')

        for RetrievalDate, URL in list(info.items()):
            schema_date = species_name + '_' + RetrievalDate
            outDit = os.path.join(pubmlst_dir, schema_date)  # compatible with windows? See if it already exists, if so, break

            if os.path.isdir(outDit):
                pickle = os.path.join(outDit, str(schema_date + '.pkl'))
                if os.path.isfile(pickle):
                    print("\tschema files already exist for %s" % (SchemaName))
                    mlst_dicts = utils.extract_variable_from_pickle(pickle)
                    SequenceDict = mlst_dicts[0]
                    for lociName, alleleSequences in list(SequenceDict.items()):
                        for sequence in alleleSequences:
                            if lociName not in list(mlst_sequences.keys()):
                                mlst_sequences[lociName] = sequence
                            else:
                                break
                    return mlst_dicts, mlst_sequences

            elif any(species_name in x for x in os.listdir(pubmlst_dir)):
                print("Older version of %s's scheme found! Deleting..." % (SchemaName))
                for directory in glob(str(pubmlst_dir + str(species_name + '_*'))):
                    utils.remove_directory(directory)
                    os.makedirs(outDit)
            else:
                os.makedirs(outDit)

            contentProfile = urllib.request.urlopen(URL[0])
            header = next(contentProfile).decode("utf8").strip().split('\t')  # skip header
            try:
                indexCC = header.index('clonal_complex') if 'clonal_complex' in header else header.index('CC')
            except:
                indexCC = len(header)
            lociOrder = header[1:indexCC]
            for row in contentProfile:
                row = row.decode("utf8").strip().split('\t')
                ST = row[0]
                alleles = ','.join(row[1:indexCC])
                STdict[alleles] = ST
            for lociName, lociURL in list(URL[1].items()):
                if lociName not in list(SequenceDict.keys()):
                    SequenceDict[lociName] = {}
                url_file = os.path.join(outDit, lociURL.rsplit('/', 1)[1])
                urllib.request.urlretrieve(lociURL, url_file)
                sequences, ignore, ignore = rematch_module.get_sequence_information(url_file, 0)
                for key in list(sequences.keys()):
                    header = re.sub("\D", "", sequences[key]['header'])
                    sequence = sequences[key]['sequence'].upper()
                    SequenceDict[lociName][sequence] = header
                    if lociName not in list(mlst_sequences.keys()):
                        mlst_sequences[lociName] = sequence
                os.remove(url_file)
            mlst_dicts = [SequenceDict, STdict, lociOrder]
            utils.save_variable_to_pickle(mlst_dicts, outDit, schema_date)
    return mlst_dicts, mlst_sequences

import sys
import os
import urllib2
import urllib
import csv
from glob import glob
import re
import functools
try:
	import xml.etree.cElementTree as ET
except ImportError:
	import xml.etree.ElementTree as ET
import utils
import rematch_module


def determine_schema(species, schema_number):
	if schema_number is not None:
		species = species + '#' + str(schema_number)
	return species.upper()


downloadPubMLST = functools.partial(utils.timer, name='Download PubMLST module')


@downloadPubMLST
def downloadPubMLSTxml(originalSpecies, schema_number, outdir):
	print '\n' + 'Searching MLST database for ' + originalSpecies

	species = determine_schema(originalSpecies, schema_number)

	xmlURL = 'http://pubmlst.org/data/dbases.xml'
	try:
		content = urllib2.urlopen(xmlURL)
		xml = content.read()
		tree = ET.fromstring(xml)
	except:
		print "Ooops! There might be a problem with the PubMLST service, try later or check if the xml is well formated at " + xmlURL
		raise

	xmlData = {}

	success = 0
	for scheme in tree.findall('species'):
		if scheme.text.strip().upper() == species or species in scheme.text.strip().upper():
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
		keys = xmlData.keys()
		keys = sorted(keys)
		print "\tWarning. More than one schema found for %s. only keeping the first one... %s" % (originalSpecies, keys[0])
		for key in keys[1:]:
			del xmlData[key]

	pubmlst_dir = os.path.join(outdir, 'pubmlst', '')
	if not os.path.isdir(pubmlst_dir):
		os.makedirs(pubmlst_dir)

	for SchemaName, info in xmlData.items():
		STdict = {}
		SequenceDict = {}
		mlst_sequences = {}

		for RetrievalDate, URL in info.items():
			schema_date = SchemaName.replace(' ', '_') + '_' + RetrievalDate
			outDit = os.path.join(pubmlst_dir, schema_date)  # compatible with windows? See if it already exists, if so, break

			if os.path.isdir(outDit):
				pickle = os.path.join(outDit, str(schema_date + '.pkl'))
				if os.path.isfile(pickle):
					print "\tschema files already exist for %s" % (SchemaName)
					mlst_dicts = utils.extractVariableFromPickle(pickle)
					SequenceDict=mlst_dicts[0]
					for lociName, alleleSequences in SequenceDict.items():
						for sequence in alleleSequences:
							if lociName not in mlst_sequences.keys():
								mlst_sequences[lociName] = sequence
							else:
								break
					return mlst_dicts, mlst_sequences

			elif any(SchemaName.replace(' ', '_') in x for x in os.listdir(pubmlst_dir)):
				print "Older version of %s's scheme found! Deleting..." % (SchemaName)
				for directory in glob(str(outDit + str(SchemaName.replace(' ', '_') + '_*'))):
					utils.removeDirectory(directory)
					os.makedirs(outDit)
			else:
				os.makedirs(outDit)

			contentProfile = urllib2.urlopen(URL[0])
			profileFile = csv.reader(contentProfile, delimiter='\t')
			header=profileFile.next()  # skip header
			try:
				indexCC=header.index('clonal_complex')
			except:
				indexCC=len(header)+1
			lociOrder=header[1:indexCC]
			for row in profileFile:
				ST = row[0]
				alleles=','.join(row[1:indexCC])
				STdict[alleles] = ST
			for lociName, lociURL in URL[1].items():
				if lociName not in SequenceDict.keys():
					SequenceDict[lociName] = {}
				url_file = os.path.join(outDit, lociURL.rsplit('/', 1)[1])
				urllib.urlretrieve(lociURL, url_file)
				sequences, headers = rematch_module.get_sequence_information(url_file, 0)
				for key in sequences.keys():
					header = re.sub("\D", "", sequences[key]['header'])
					sequence = sequences[key]['sequence'].upper()
					SequenceDict[lociName][sequence] = header
					if lociName not in mlst_sequences.keys():
						mlst_sequences[lociName] = sequence
				os.remove(url_file)
			mlst_dicts = [SequenceDict, STdict, lociOrder]
			utils.saveVariableToPickle(mlst_dicts, outDit, schema_date)
	return mlst_dicts, mlst_sequences

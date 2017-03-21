import sys
import os
import urllib2
import urllib
import csv
from glob import glob
import re
try:
	import xml.etree.cElementTree as ET
except ImportError:
	import xml.etree.ElementTree as ET
import utils
import rematch_module


def downloadPubMLSTxml(originalSpecies, outdir):
	xmlURL = 'http://pubmlst.org/data/dbases.xml'
	print '\n' + 'Searching MLST database for ' + originalSpecies

	species = originalSpecies.upper()

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

		for RetrievalDate, URL in info.items():
			outDit = os.path.join(pubmlst_dir, str(SchemaName.replace(' ', '_') + '_' + RetrievalDate))  # compatible with windows? See if it already exists, if so, break

			if os.path.isdir(outDit):
				pickle = os.path.join(outDit, str(SchemaName.replace(' ', '_') + '_' + RetrievalDate + '.pkl'))
				if os.path.isfile(pickle):
					print "\tschema files already exist for %s" % (SchemaName)
					toSave = utils.extractVariableFromPickle(pickle)
					break
				else:
					print 'MPM: 1'
			elif any(SchemaName.replace(' ', '_') in x for x in os.listdir(pubmlst_dir)):
				print "Older version of %s's scheme found! Deleting..." % (SchemaName)
				for directory in glob(str(outDit + str(SchemaName.replace(' ', '_') + '_*'))):
					utils.removeDirectory(directory)
					os.makedirs(outDit)
			else:
				os.makedirs(outDit)

			contentProfile = urllib2.urlopen(URL[0])
			profileFile = csv.reader(contentProfile, delimiter='\t')
			profileFile.next()  # skip header
			for row in profileFile:
				ST = row[0]
				alleles = row[1:-1]
				STdict[','.join(alleles)] = ST
			for lociName, lociURL in URL[1].items():
				if lociName not in SequenceDict.keys():
					SequenceDict[lociName] = {}
				urllib.urlretrieve(lociURL, os.path.join(outDit, lociURL.rsplit('/', 1)[1]))
				sequences, headers=rematch_module.get_sequence_information(outDit+'/'+os.path.basename(lociURL),0)
				for key in sequences.keys():
					header=re.sub("\D", "", sequences[key]['header'])
					sequence=sequences[key]['sequence'].upper()
					SequenceDict[lociName][sequence]=header
				os.remove(outDit+'/'+os.path.basename(lociURL))
			toSave=[SequenceDict,STdict]
			utils.saveVariableToPickle(toSave,outDit,SchemaName.replace(' ','_')+'_'+RetrievalDate)
	'''
	for k,v in STdict.items():
		print v + '->' + k

	for k,v in SequenceDict.items():
		print k
		for key,value in v.items():
			print value
			print key

	if not ToSkip:
		toSave=[SequenceDict,STdict]
		utils.saveVariableToPickle(toSave,outDit,SchemaName.replace(' ','_')+'_'+RetrievalDate)
		return True
	else:
		return False

	'''
	#print SequenceDict, STdict
	#print len(out)
	return toSave

def main():
	outdir = '.'
	'''agalactiae=downloadPubMLSTxml('Streptococcus agalactiae')

	pyogenes=downloadPubMLSTxml('Streptococcus pyogenes')

	pneumo=downloadPubMLSTxml('Streptococcus pneumoniae')

	aspergillus=downloadPubMLSTxml('Aspergillus fumigatus')
	'''
	salmonella=downloadPubMLSTxml('Salmonella enterica', outdir)

	yersinia=downloadPubMLSTxml('Yersinia pseudotuberculosis', outdir)

	campy=downloadPubMLSTxml('Campylobacter jejuni', outdir)

	acineto, lala=downloadPubMLSTxml('Acinetobacter baumannii', outdir)
	#print acineto
	#print lala

	downloadPubMLSTxml('lala', outdir) #TODO

if __name__ == "__main__":
	main()

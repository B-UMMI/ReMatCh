import argparse, sys, os
import urllib2
import urllib
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import time
import csv
import utils
import rematch_module
import utils
from  __builtin__ import any
from glob import glob
import shutil

def downloadPubMLSTxml(species, xmlURL='http://pubmlst.org/data/dbases.xml'):

	print '\n' + 'Searching MLST database for ' + species

	try:
		content = urllib2.urlopen(xmlURL)
		xml = content.read()
		tree = ET.fromstring(xml)
		#print tree #Element Root
	except:
		print "Ooops!There might be a problem with the PubMLST service, try later or check if the xml is well formated at " + xmlURL
		raise

	xmlData={}
	for scheme in tree.findall('species'):
		try:
			if scheme.text.strip() == species:
				xmlData[scheme.text.strip()]={}
				for info in scheme: #mlst
					for database in info: #database
						for retrievedDate in database.findall('retrieved'):
							retrieved = retrievedDate.text
							xmlData[scheme.text.strip()][retrieved]=[]
						for profile in database.findall('profiles'):
								profileURl=profile.find('url').text
								xmlData[scheme.text.strip()][retrieved].append(profileURl)
						for lociScheme in database.findall('loci'):
							loci={}
							for locus in lociScheme:
								locusID=locus.text
								for locusInfo in locus:
									locusUrl=locusInfo.text
									loci[locusID.strip()]=locusUrl
								xmlData[scheme.text.strip()][retrieved].append(loci)
		except:
			print "No schema found for %s" % (species) #TODO - not working
			raise

	#uncomment for verbose
	'''
	print '\n'
	for key,value in xmlData.items():
		print 'species: ' + key
		for k,v in value.items():
			print 'Date of Retrieval: ' + k
			print 'Profile URL: ' + v[0]
			for loci, url in v[1].items():
				print "Loci name: " + loci
				print "Loci URL: " + url
	'''
	return xmlData

def getSequences(URLdictionary):

	ToSkip=False
	STdict={}
	SequenceDict={}

	if not os.path.isdir('./pubmlst'):
		os.makedirs('./pubmlst')

	for SchemaName, info in URLdictionary.items():
		#print SchemaName
		for RetrievalDate, URL in info.items():
			#print RetrievalDate
			outDit='./pubmlst'+'/'+SchemaName.replace(' ','_')+'_'+RetrievalDate #compatible with windows? See if it already exists, if so, break
			if os.path.isdir(outDit):
				print "\tschema files already exist for %s" % (SchemaName)
				ToSkip=True
				break
			elif any(SchemaName.replace(' ','_') in x for x in os.listdir('./pubmlst/')):
				print "Older version of %s's scheme found! Deleting..." % (SchemaName)
				for directory in glob(str('./pubmlst'+'/'+SchemaName.replace(' ','_')+'_*')):
					shutil.rmtree(directory)
					os.makedirs(outDit)
			else:
				os.makedirs(outDit)
			contentProfile=urllib2.urlopen(URL[0])
			profileFile =  csv.reader(contentProfile, delimiter='\t')
			profileFile.next() #skip header
			for row in profileFile:
				ST=row[0]
				alleles=row[1:-1]
				STdict[','.join(alleles)]=ST 
			for lociName, lociURL in URL[1].items():
				if lociName not in SequenceDict.keys():
					SequenceDict[lociName]={}
				contentLoci=urllib.urlretrieve(lociURL, outDit+'/'+os.path.basename(lociURL))
				sequences, headers=rematch_module.get_sequence_information(outDit+'/'+os.path.basename(lociURL),0)
				for key in sequences.keys():
					header=sequences[key]['header'].split('_')[-1]
					sequence=sequences[key]['sequence']
					SequenceDict[lociName][sequence]=header
				os.remove(outDit+'/'+os.path.basename(lociURL))
	'''
	for k,v in STdict.items():
		print v + '->' + k

	for k,v in SequenceDict.items():
		print k
		for key,value in v.items():
			print value
			print key
	'''
	if not ToSkip:
		toSave=[SequenceDict,STdict]
		utils.saveVariableToPickle(toSave,outDit,SchemaName.replace(' ','_')+'_'+RetrievalDate)







def main():
	agalactiae=downloadPubMLSTxml('Streptococcus agalactiae')
	getSequences(agalactiae)

	pyogenes=downloadPubMLSTxml('Streptococcus pyogenes')
	getSequences(pyogenes)

	pneumo=downloadPubMLSTxml('Streptococcus pneumoniae')
	getSequences(pneumo)

	aspergillus=downloadPubMLSTxml('Aspergillus fumigatus')
	getSequences(aspergillus)

	salmonella=downloadPubMLSTxml('Salmonella enterica')
	getSequences(salmonella)

	campy=downloadPubMLSTxml('Campylobacter jejuni')
	getSequences(campy)
	
	downloadPubMLSTxml('lala') #TODO

if __name__ == "__main__":
    main()


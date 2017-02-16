#!/usr/bin/env python

import argparse, sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import ntpath

def parseID(filename):
	#get wanted feature IDs
	gffIDs=[]
	with open(filename, 'r') as in_handle:
		for line in in_handle:
			line=line.strip()
			gffIDs.append(line)
	return gffIDs

def retrieveSeq_File(fastaFile, coordFile, extraSeq, filename, outputDir):
	#Parsing the sequence file, using the provided txt file containing the contig ID and positions to retrieve sequences.
	handle = open(fastaFile, "rU")
	records_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	handle.close()

	Seq2Get={}
	with open(coordFile, 'r') as sequeces2get:
		for line in sequeces2get:
			line=line.split(',')
			coords=(int(line[-2]), int(line[-1]))
			contigID=line[0]
			if contigID in Seq2Get.keys():
				Seq2Get[contigID].append(coords)
			else:
				Seq2Get[contigID]=[coords]

	with open(outputDir+'/'+filename+'.fasta','w') as output_handle:
		fails=0
		successes=0
		records=[]
		for contig, listCoords in Seq2Get.items():
			contigSeq=records_dict[contig].seq
			for coord in listCoords:
				coord1=coord[0]-extraSeq
				coord2=coord[1]+extraSeq
				if coord1 < 0 or coord2 > len(contigSeq):
					fail_log=open(outputDir+'/'+filename+'_fails.txt', 'a')
					fail_log.write(contig + ',' + str(coord[0])+','+ str(coord[1])+'\n')
					fail_log.close()
					fails+=1
				else:
					geneseq=str(contigSeq[coord1:coord2])
					record = SeqRecord(Seq(geneseq), id=str(str(contig)+'#'+str(coord1)+'_'+str(coord2)), description='')
					records.append(record)
					successes+=1
		SeqIO.write(records, output_handle, "fasta")

	print 'Retrived %s features successfully from %s with %s bp as extra sequence.' % (str(successes), filename, str(extraSeq))
	if fails>0:
		print '%s featrued failed to retrieve. Check %s_fails.txt file.' % (str(fails), filename)

def retrieveSeq(fastaFile, gffFeatures, extraSeq, filename, outputDir):
	#parsing the sequence file into a SeqIO dictionary. one contig per entry
	handle = open(fastaFile, "rU")
	records_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	handle.close()

	with open(outputDir+'/'+filename+'.fasta','w') as output_handle:
		fails=0
		successes=0
		records=[]
		for locus, location in gffFeatures.items():
			#print locus
			contigSeq=records_dict[location[0]].seq
			coord1=location[1]-extraSeq
			coord2=location[2]+extraSeq
			if coord1 < 0 or coord2 > len(contigSeq):
				fail_log=open(outputDir+'/'+filename+'_fails.txt', 'a')
				fail_log.write(locus+'\n')
				fail_log.close()
				fails+=1
			else:
				geneseq=str(contigSeq[coord1:coord2])
				if location[3] == '-':
					seq = Seq(geneseq)
					geneseq =str(seq.reverse_complement())
				record = SeqRecord(Seq(geneseq), id=str(locus+'-'+str(location[0])+'#'+str(location[1])+'_'+str(location[2])), description='')
				records.append(record)
				successes+=1
		SeqIO.write(records, output_handle, "fasta")
	print 'Retrived %s features successfully from %s with %s bp as extra sequence.' % (str(successes), filename, str(extraSeq))
	if fails>0:
		print '%s featrued failed to retrieve. Check %s_fails.txt file.' % (str(fails), filename)

def parseFeatures(tempGFF):
	#parsing the feature file into a dictionary
	gffFeatures={}

	with open(tempGFF, 'r') as temp_genes:
		for line in temp_genes:
			line=line.split('\t')
			if "CDS" in line[2]:
				ID=line[-1].split(';')
				locusID=str(ID[0].split('=')[1])
				contig=line[0]
				begining=int(line[3])-1 #to get the full sequence
				end=int(line[4])
				strand=line[6]
				location=[contig, begining, end, strand]
				gffFeatures[locusID]=location
	return gffFeatures

def gffParser(gffFile, extraSeq=0, outputDir='.', keepTemporaryFiles=False, IDs=None, coordFile=None):

	filename=ntpath.basename(gffFile).replace('.gff', '')

	#cleaning temp files if they exist
	if os.path.isfile(outputDir+'/'+filename+'_features.gff'):
		os.remove(outputDir+'/'+filename+'_features.gff')
	if os.path.isfile(outputDir+'/'+filename+'_sequence.fasta'):
		os.remove(outputDir+'/'+filename+'_sequence.fasta')

	#cleaning fails file if it exists
	if os.path.isfile(outputDir+'/'+filename+'_fails.txt'):
		os.remove(outputDir+'/'+filename+'_fails.txt')

	if coordFile is None:

		if IDs is not None:
			selectIDs=parseID(IDs)
		else:
			selectIDs=None
	
		#separating the gff into 2 different files: one with the features and another with the conting sequences
		with open(gffFile, 'r') as in_handle, open(outputDir+'/'+filename+'_features.gff', 'a') as temp_genes, open(outputDir+'/'+filename+'_sequence.fasta', 'a') as temp_contigs:
			for line in in_handle: 
				if not line.startswith('##'):
					if '\t' in line:
						if selectIDs is not None:
							items=line.split('\t')
							ID=items[-1].split(';')[0]
							ID=ID.split('=')[1]
							if ID in selectIDs:
								temp_genes.write(line)
						else:
							temp_genes.write(line)
					else:
						temp_contigs.write(line)

		gffFiles=parseFeatures(outputDir+'/'+filename+'_features.gff')

		retrieveSeq(outputDir+'/'+filename+'_sequence.fasta', gffFiles, extraSeq, filename, outputDir)
	
	else:
		with open(gffFile, 'r') as in_handle, open(outputDir+'/'+filename+'_sequence.fasta', 'a') as temp_contigs:
			for line in in_handle:
				if not line.startswith('##'):
					if '\t' in line:
						pass
					else:
						temp_contigs.write(line)

		retrieveSeq_File(outputDir+'/'+filename+'_sequence.fasta', coordFile, extraSeq, filename, outputDir)

	#removing temp files
	if not keepTemporaryFiles:
		try:
			os.remove(outputDir+'/'+filename+'_features.gff')
		except:
			pass
		os.remove(outputDir+'/'+filename+'_sequence.fasta')

def main():

	version='1.0.0'

	parser = argparse.ArgumentParser(prog='gffParser.py', description='GFF3 parser for feature sequence retrival.', epilog='by C I Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	parser.add_argument('-i', '--input', help='GFF3 file to parse, containing both sequences and annotations (like the one obtained from PROKKA).', type=argparse.FileType('r'), required=True)
	parser.add_argument('-x', '--extraSeq', help='Extra sequence to retrieve per feature in gff.', default=0, type=int, required=False)
	parser.add_argument('-k','--keepTemporaryFiles', help='Keep temporary gff(without sequence) and fasta files.', action='store_true')
	parser.add_argument('-o', '--outputDir', help='Path to where the output is to be saved.', default='.', required=False)
	

	parser_optional_selected_regions_exclusive = parser.add_mutually_exclusive_group()
	parser_optional_selected_regions_exclusive.add_argument('-s', '--select', help='txt file with the IDs of interest, one per line', type=argparse.FileType('r'), required=False)
	parser_optional_selected_regions_exclusive.add_argument('-f', '--fromFile', help='Sequence coordinates to be retrieved. Requires contig ID and coords (contig,strart,end) in a csv file. One per line.', type=argparse.FileType('r'), required=False)

	args = parser.parse_args()
	
	gffParser(os.path.abspath(args.input.name), args.extraSeq, os.path.abspath(args.outputDir), args.keepTemporaryFiles, os.path.abspath(args.select.name) if args.select is not None else None, os.path.abspath(args.fromFile.name) if args.fromFile is not None else None)


if __name__ == "__main__":
    main()

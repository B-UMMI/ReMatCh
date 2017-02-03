#!/usr/bin/env python

import argparse, sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import ntpath


def parseID(filename):
	gffIDs=[]
	with open(filename, 'r') as in_handle:
		for line in in_handle:
			line=line.strip()
			gffIDs.append(line)
	return gffIDs



def retrieveSeq(fastaFile, gffFeatures, extraSeq, filename, outputDir):
	#parsing the sequence file into a SeqIO dictionary. one contig per entry
	handle = open(fastaFile, "rU")
	records_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	handle.close()

	with open(outputDir+'/'+filename+'.fasta','w') as output_handle:
		lala=0
		records=[]
		for locus, location in gffFeatures.items():
			#print locus
			contigSeq=records_dict[location[0]].seq
			coord1=location[1]-extraSeq
			coord2=location[2]+extraSeq
			if coord1 < 0 or coord2 > len(contigSeq):
				print locus
				lala+=1
			else:
				geneseq=str(contigSeq[coord1:coord2])
				if location[3] == '-':
					seq = Seq(geneseq)
					geneseq =str(seq.reverse_complement())
				record = SeqRecord(Seq(geneseq), id=str(locus+':'+str(location[0])+'_'+str(location[1])+'_'+str(location[2])), description='')
				records.append(record)
		SeqIO.write(records, output_handle, "fasta")
		print lala
			#output.write('>'+locus+'_'+str(location[0])+'_'+str(location[1])+'_'+str(location[2])+'\n')
			#output.write(geneseq)

def parseFeatures(tempGFF):
	#parsing the feature file into a dictionary
	gffFeatures={}

	with open(tempGFF, 'r') as temp_genes:
		for line in temp_genes:
			line=line.split('\t')
			#print line
			if "CDS" in line[2]:
				ID=line[-1].split(';')
				locusID=str(ID[0].split('=')[1])
				#print locusID
				contig=line[0]
				begining=int(line[3])-1 #to get the full sequence
				end=int(line[4])
				strand=line[6]
				location=[contig, begining, end, strand]
				gffFeatures[locusID]=location
	return gffFeatures

def gffParser(gffFile, extraSeq, outputDir, keepTemporaryFiles, selectIDs):

	filename=ntpath.basename(gffFile).replace('.gff', '')
	#print filename

	#cleaning temp files if they exist
	if os.path.isfile(outputDir+'/'+filename+'_features.gff'):
		os.remove(outputDir+'/'+filename+'_features.gff')
	if os.path.isfile(outputDir+'/'+filename+'_sequence.fasta'):
		os.remove(outputDir+'/'+filename+'_sequence.fasta')
	
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
	

	#removing temp files
	if not keepTemporaryFiles:
		os.remove(outputDir+'/'+filename+'_features.gff')
		os.remove(outputDir+'/'+filename+'_sequence.fasta')




def main():

	version='1.4.5'

	parser = argparse.ArgumentParser(description='GFF3 parser for feature sequence retrival, containing both sequences and annotations.', epilog='by C I Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-i', '--input', help='GFF3 file to parse, containing both sequences and annotations (like the one obtained from PROKKA).', required=True)#,type=type=argparse.FileType('r'), required=True)
	parser.add_argument('-x', '--extra_seq', help='Extra sequence to retrieve per feature in gff.', default=0, type=int, required=False)
	parser.add_argument('-k','--keepTemporaryFiles', help='Keep temporary gff(without sequence) and fasta files.', default=False, action='store_true', required=False)
	parser.add_argument('-o', '--outputDir', help='Path to where the output is to be saved.', default='.')
	parser.add_argument('-s', '--select', help='txt file with the IDs of interest, one per line', default=None)
	#parser.add-argument('--verbose')
	parser.add_argument('--version', help='Display version, and exit.', default=False, action='store_true')

	args = parser.parse_args()

	#version
	if args.version:
		print sys.stdout, "Current version: %s" %(version)
		sys.exit(0)

	#check if something is missing
	if args.input is None:
		parser.print_usage()
		print "error: argument -i/--input is required"
		sys.exit(1)

	#START

	if args.select is not None:
		IDs=parseID(args.select)
	else:
		IDs=None

	print 'parsing gff file...'
	gffParser(args.input, args.extra_seq, args.outputDir, args.keepTemporaryFiles, IDs)

	print "Finished"
	sys.exit(0)


if __name__ == "__main__":
    main()
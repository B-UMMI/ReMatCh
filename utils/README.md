ReMatCh
=======
*Reads mapping against target sequences, checking mapping and consensus sequences production*  

<https://github.com/B-UMMI/ReMatCh>


gffParser
---------

Parser for GFF3 files, as the ones obtained by [PROKKA] (https://github.com/tseemann/prokka). This files require to have both the features and sequence. It will retrieve the CDS sequences in the GFF file, allowing these to be extended by the number of nucleotides specifiend in `--extraSeq`. A selection of CDS of interest to be parsed can also be obtained by providing `--select` with a txt file of the IDs of interest, one per line. As an alternative, wanted sequences can be obtained from the GFF file from a txt file containing the coontig ID, start and end position (one per line) of the sequences of interest, using the `-fromFile` option. `-extraSeq` can also be obtain through this method.

**Dependencies**
- Python (2.7.x)
- [Biopython] (http://biopython.org/) (1.68 or similar)

**Usage**

    gffParser.py 	[-h]
    				-i INPUT [-x EXTRASEQ] [-k] [-o OUTPUTDIR]
	                [-s SELECT] [-f FROMFILE] [--version]

	GFF3 parser for feature sequence retrival, containing both sequences and
	annotations.

	optional arguments:
  		-h, --help          Show this help message and exit
  		-i INPUT, --input INPUT
                        	GFF3 file to parse, containing both sequences and
                        	annotations (like the one obtained from PROKKA).
  		-x EXTRASEQ, --extraSeq EXTRASEQ
                        	Extra sequence to retrieve per feature in gff.
  		-k, --keepTemporaryFiles
                        	Keep temporary gff(without sequence) and fasta files.
  		-o OUTPUTDIR, --outputDir OUTPUTDIR
                        	Path to where the output is to be saved.
  		-s SELECT, --select SELECT
                        	txt file with the IDs of interest, one per line
  		-f FROMFILE, --fromFile FROMFILE
                        	Sequence coordinates to be retrieved. Requires contig ID
                        	and coords (contig,strart,end) in a csv file, one per line.
  		--version           Display version, and exit.

**Output**

*<filename>.fasta*  
Multi-fasta file with the retrieved sequences.
Headers will contain the feature ID, followed by '=', and the position of that feature in the sequence, starting with the original sequence ID,a '# and' the start and end coordinates separated with '_' (>featureID=contig#start_end).
If the `--fromFile` option is used, there's no feature ID, so the header will only contain it's position in the original sequence, followed by the start and end coordinates separated with '_' (>contig#start_end).

*<filename>.txt*  
Feature ID of the sequences that failed to be retireved, due to the start position or end position being outside of the sequence where the feature is (due to the `--extraSeq` option).


restart_rematch
---------------

Restart a ReMatCh run abruptly terminated

**Usage**

  usage: restart_rematch.py [-h] [--version] -i
                            /path/to/initial/workdir/directory/
                            [-w /path/to/workdir/directory/] [-j N]
                            [--runFailedSamples]

  Restart a ReMatCh run abruptly terminated

  optional arguments:
    -h, --help            show this help message and exit
    --version             Version information

  Required options:
    -i /path/to/initial/workdir/directory/, --initialWorkdir /path/to/initial/workdir/directory/
                          Path to the directory where ReMatCh was running
                          (default: None)

  General facultative options:
    -w /path/to/workdir/directory/, --workdir /path/to/workdir/directory/
                          Path to the directory where ReMatCh will run again
                          (default: .)
    -j N, --threads N     Number of threads to use instead of the ones set in
                          initial ReMatCh run (default: None)
    --runFailedSamples    Will run ReMatCh for those samples missing, as well as
                          for samples that did not run successfully in initial
                          ReMatCh run (default: False)
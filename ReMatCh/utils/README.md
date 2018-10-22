ReMatCh
=======
*Reads mapping against target sequences, checking mapping and consensus sequences production*  

<https://github.com/B-UMMI/ReMatCh>

Table of Contents
--

[Combine alignment consensus](#combine-alignment-consensus)
[Convert Ns to gaps](#convert-ns-to-gaps)
[gffParser](#gffparser)
[Restart ReMatCh](#restart-rematch)
[Strip Alignment](#strip-alignment)


## Combine Alignment Consensus

Combine the alignment consensus sequences from ReMatCh first run by reference sequences into single files.

**Dependencies**
- Python v3

**Usage**

    usage: combine_alignment_consensus.py [-h] [--version]
                      -w                 /path/to/rematch/working/directory/
                      [-o /path/to/output/directory/]

    Combine the alignment consensus sequences from ReMatCh first run by reference sequences into single files

    optional arguments:
                    -h, --help show this help message and exit
                    --version Version information

    Required options:
                    -w /path/to/rematch/working/directory/
                    --workdir /path/to/rematch/working/directory/ Path to the directory where ReMatCh was running (default: None)

    General facultative options:
                    -o --outdir /path/to/output/directory/  Path to the directory where the combined sequence files will stored (default: .)



## Convert Ns to Gaps


Convert the Ns into gaps in a fasta file.

**Dependencies**
- Python (2.7.x)

**Usage**

    usage: convert_Ns_to_gaps.py [-h] [--version]
                    -i /path/to/input/file.fasta
                    -o                     /path/to/converted/output/file.fasta

    Convert the Ns into gaps

    optional arguments:
                    -h, --help show this help message and exit
                    --version Version information

    Required options:
                    -i --infile /path/to/input/file.fasta Path to the fasta file (default: None)
                    -o --outfile  /path/to/converted/output/file.fasta                        Converted output fasta file (default:                        converted_Ns_to_gaps.fasta)



## gffParser


Parser for GFF3 files, as the ones obtained by [PROKKA](https://github.com/tseemann/prokka). This files require to have both the features and sequence. It will retrieve the CDS sequences in the GFF file, allowing these to be extended by the number of nucleotides specifiend in `--extraSeq`. A selection of CDS of interest to be parsed can also be obtained by providing `--select` with a txt file of the IDs of interest, one per line. As an alternative, wanted sequences can be obtained from the GFF file from a txt file containing the coontig ID, start and end position (one per line) of the sequences of interest, using the `-fromFile` option. `-extraSeq` can also be obtain through this method.

**Dependencies**
- Python (2.7.x)
- [Biopython](http://biopython.org/) (1.68 or similar)

**Usage**

    usage: gffParser.py 	[-h]
                    -i INPUT [-x EXTRASEQ] [-k] [-o OUTPUTDIR]
	                  [-s SELECT] [-f FROMFILE] [--version]

	  GFF3 parser for feature sequence retrival, containing both sequences and annotations.

	optional arguments:
                    -h, --help    Show this help message and exit
                    -i --input INPUT
                                  GFF3 file to parse, containing both sequences and annotations (like the one obtained from PROKKA).
                    -x --extraSeq EXTRASEQ
                                  Extra sequence to retrieve per feature in gff.
                    -k, --keepTemporaryFiles
                                  Keep temporary gff(without sequence) and fasta files.
                    -o --outputDir OUTPUTDIR
                                  Path to where the output is to be saved.
                    -s --select SELECT
                                  txt file with the IDs of interest, one per line
                    -f --fromFile FROMFILE
                                  Sequence coordinates to be retrieved. Requires contig ID and coords (contig,strart,end) in a csv file, one per line.
                    --version     Display version, and exit.

**Output**

*<filename>.fasta*  
Multi-fasta file with the retrieved sequences.
Headers will contain the feature ID, followed by '=', and the position of that feature in the sequence, starting with the original sequence ID,a '# and' the start and end coordinates separated with '_' (>featureID=contig#start_end).
If the `--fromFile` option is used, there's no feature ID, so the header will only contain it's position in the original sequence, followed by the start and end coordinates separated with '_' (>contig#start_end).

*<filename>.txt*  
Feature ID of the sequences that failed to be retireved, due to the start position or end position being outside of the sequence where the feature is (due to the `--extraSeq` option).



## Restart ReMatCh


Restart a ReMatCh run abruptly terminated

**Dependencies**
- Python (2.7.x)

**Usage**

    usage: restart_rematch.py [-h] [--version] -i
                            /path/to/initial/workdir/directory/
                            [-w /path/to/workdir/directory/] [-j N]
                            [--runFailedSamples]

    Restart a ReMatCh run abruptly terminated

    optional arguments:
                    -h, --help    show this help message and exit
                    --version     Version information

    Required options:
                    -i /path/to/initial/workdir/directory/, --initialWorkdir /path/to/initial/workdir/directory/
                                  Path to the directory where ReMatCh was running (default: None)

  General facultative options:
                    -w, --workdir /path/to/workdir/directory/
                                  Path to the directory where ReMatCh will run again (default: .)
                    -j N, --threads N
                                  Number of threads to use instead of the ones set in initial ReMatCh run (default: None)
                    --runFailedSamples
                                  Will run ReMatCh for those samples missing, as well as for samples that did not run successfully in initial ReMatCh run (default: False)



## Strip Alignment


Strip alignment positions containing gaps,
missing data and invariable positions.

**Dependencies**
- Python (2.7.x)
- [Biopython](http://biopython.org/) (1.68 or similar)

**Usage**

    usage: strip_alignment.py [-h] [--version]
                    -i                      /path/to/aligned/input/file.fasta -o /path/to/stripped/output/file.fasta [--notGAPs]
                    [--notMissing] [--notInvariable]

    Strip alignment positions containing gaps, missing data and invariable positions

    optional arguments:
                    -h, --help    show this help message and exit
                    --version     Version information

    Required options:
                    -i, --infile /path/to/aligned/input/file.fasta
                                  Path to the aligned fasta file (default: None)
                    -o, --outfile /path/to/stripped/output/file.fasta
                                  Stripped output fasta file (default: alignment_stripped.fasta)

    General facultative options:
                    --notGAPs     Not strip positions with GAPs (default: False)
                    --notMissing  Not strip positions with missing data (default: False)
                    --notInvariable
                                  Not strip invariable sites (default: False)

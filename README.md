ReMatCh
=======
*Reads mapping against target sequences, checking mapping and consensus sequences production*  

<https://github.com/B-UMMI/ReMatCh>



Requirements
------------

 - Fasta file with reference sequences



Dependencies
------------
**Optional**

Required to download sequence data from ENA database:
 - *Aspera Connect 2* >= v3.6.1
 - *gzip* >= v1.6 (normally found in Linux OS) (important for bam/cram conversion into fastq)

Required to run ReMatch analysis
 - *Bowtie2* >= v2.2.9
 - *Samtools* = v1.3.1
 - *Bcftools* = v1.3.1  
(these three executables are provided, but user's own executables can be used by providing `--doNotUseProvidedSoftware` option)



Installation
------------
    git clone https://github.com/B-UMMI/ReMatCh.git



Usage
-----
    usage: rematch.py [-h] [--version]
                      (-r /path/to/reference_sequence.fasta | --mlstReference)
                      [-w /path/to/workdir/directory/] [-j N]
                      [--mlst "Streptococcus agalactiae"]
                      [--doNotUseProvidedSoftware] [--extraSeq N]
                      [--minCovPresence N] [--minCovCall N]
                      [--minFrequencyDominantAllele 0.6] [--minGeneCoverage N]
                      [--minGeneIdentity N] [--doubleRun] [--notWriteConsensus]
                      [--bowtieOPT] [--debug]
                      [--mlstSchemaNumber N] [--mlstConsensus noMatter]
                      [--mlstRun first]
                      [-a /path/to/asperaweb_id_dsa.openssh] [-k]
                      [--downloadLibrariesType PAIRED]
                      [--downloadInstrumentPlatform ILLUMINA] [--downloadCramBam]
                      [-l /path/to/list_IDs.txt | -t "Streptococcus agalactiae"]

    Reads mapping against target sequences, checking mapping and consensus
    sequences production

    optional arguments:
      -h, --help            show this help message and exit
      --version             Version information
      -l /path/to/list_IDs.txt, --listIDs /path/to/list_IDs.txt
                      Path to list containing the IDs to be downloaded (one
                      per line) (default: None)
      -t "Streptococcus agalactiae", --taxon "Streptococcus agalactiae"
                      Taxon name for which ReMatCh will download fastq files
                      (default: None)

    General facultative options:
      -r /path/to/reference_sequence.fasta, --reference /path/to/reference_sequence.fasta
                            Fasta file containing reference sequences (default:
                            None)
      -w /path/to/workdir/directory/, --workdir /path/to/workdir/directory/
                            Path to the directory where ReMatCh will run and
                            produce the outputs with reads (ended with
                            fastq.gz/fq.gz and, in case of PE data, pair-end
                            direction coded as _R1_001 / _R2_001 or _1 / _2)
                            already present (organized in sample folders) or
                            to be downloaded (default: .)
      -j N, --threads N     Number of threads to use (default: 1)
      --mlst "Streptococcus agalactiae"
                            Species name (same as in PubMLST) to be used in MLST
                            determination (default: None)
      --doNotUseProvidedSoftware
                            Tells ReMatCh to not use Bowtie2, Samtools and
                            Bcftools that are provided with it (default: False)

    ReMatCh module facultative options:
      --extraSeq N          Sequence length added to both ends of target sequences
                            (usefull to improve reads mapping to the target one)
                            that will be trimmed in ReMatCh outputs (default: 0)
      --minCovPresence N    Reference position minimum coverage depth to consider
                            the position to be present in the sample (default: 5)
      --minCovCall N        Reference position minimum coverage depth to perform a
                            base call. Lower coverage will be coded as N (default:
                            10)
      --minFrequencyDominantAllele 0.6
                            Minimum relative frequency of the dominant allele
                            coverage depth (value between [0, 1]). Positions with
                            lower values will be considered as having multiple
                            alleles (and will be coded as N) (default: 0.6)
      --minGeneCoverage N   Minimum percentage of target reference gene sequence covered
                            by --minCovPresence to consider a gene to be present
                            (value between [0, 100]) (default: 80)
      --minGeneIdentity N   Minimum percentage of identity of reference gene sequence
                            covered by --minCovCall to consider a gene to be present
                            (value between [0, 100]). One INDEL will be considered
                            as one difference (default: 70)
      --doubleRun           Tells ReMatCh to run a second time using as reference the
                            noMatter consensus sequence produced in the first run.
                            This will improve consensus sequence determination for
                            sequences with high percentage of target reference gene
                            sequence covered (default: False)
      --notWriteConsensus   Do not write consensus sequences (default: False)
      --bowtieOPT "--no-mixed"
                            Extra Bowtie2 options (default: None)
      --debug               DeBug Mode: do not remove temporary files (default: False)
      --mlstReference       If the curated scheme for MLST alleles is available, tells
                            ReMatCh to use these as reference (force Bowtie2 to run
                            with very-sensitive-local parameters, and sets --extraSeq
                            to 200), otherwise ReMatCh uses the first alleles of each
                            MLST gene fragment in PubMLST as reference sequences (force
                            Bowtie2 to run with very-sensitive-local parameters, and
                            sets --extraSeq to 0)

    MLST facultative options:
      --mlstSchemaNumber N  Number of the species PubMLST schema to be used in
                            case of multiple schemes available (by default will
                            use the first schema) (default: None)
      --mlstConsensus noMatter
                            Consensus sequence to be used in MLST determination
                            (default: noMatter)
      --mlstRun first       ReMatCh run outputs to be used in MLST determination
                            (default: all)

    Download facultative options:
      -a /path/to/asperaweb_id_dsa.openssh, --asperaKey /path/to/asperaweb_id_dsa.openssh
                            Tells ReMatCh to download fastq files from ENA using
                            Aspera Connect. With this option, the path to Private-
                            key file asperaweb_id_dsa.openssh must be provided
                            (normaly found in
                            ~/.aspera/connect/etc/asperaweb_id_dsa.openssh).
                            (default: None)
      -k, --keepDownloadedFastq
                            Tells ReMatCh to keep the fastq files downloaded
                            (default: False)
      --downloadLibrariesType PAIRED
                            Tells ReMatCh to download files with specific library
                            layout (default: BOTH)
      --downloadInstrumentPlatform ILLUMINA
                            Tells ReMatCh to download files with specific library
                            layout (default: ILLUMINA)
      --downloadCramBam     Tells ReMatCh to also download cram/bam files and
                            convert them to fastq files (default: False)


    Download list options (one of the following):
      -l /path/to/list_IDs.txt, --listIDs /path/to/list_IDs.txt
                            Path to list containing the IDs to be downloaded (one
                            per line) (default: None)
      -t "Streptococcus agalactiae", --taxon "Streptococcus agalactiae"
                            Taxon name for which ReMatCh will download fastq files
                            (default: None)



**Running ReMatCh in local samples**  
To run ReMatCh in local fastq files, please organize those files in sample folders.  
It is advisable to use copied fastq files or symbolic links to the original files.  
Then provide the directory containing sample folders to `--workdir`. ReMatCh will store the output files there.  
E.g.:  
```
  workir/
    sample_1/
      fastq_file_a_1.fq.gz
      fastq_file_a_2.fq.gz
    sample_2/
      fastq_file_b_R1_001.fastq.gz
      fastq_file_b_R2_001.fastq.gz
```

**Running ReMatCh in specific ENA sequencing data**  
To run ReMatCh in a specific set of ENA IDs, provide a file to `--listIDs` containing a list of ENA IDs that will be downloaded.  
The IDs can be Sample Accession numbers or Run Accession numbers (for example).  
In case of IDs containing more than one Run Accession number (like Study accession numbers), only one of them will be downloaded and the remaining will be stored in *sample_report.*.tab* file under extra_run_accession column in a comma separated style.  
ReMatCh will store the output files in the `--workdir`.  

**Running ReMatCh in ENA sequencing data of a given taxon**  
To run ReMatCh in all ENA data of a given taxon, provide the taxon name to `--taxon`.  
The ENA Run Accession numbers for the given taxon will be stored in IDs_list.seqFromWebTaxon.tab file.  
The column content will be: 1) Run Accession numbers, 2) Sequencing instrument models, 3) (secondary) Study Accession numbers, 4) library types, 5) library layouts.  
The first line of *IDs_list.seqFromWebTaxon.tab* will contain the date of accession.  

**Running ReMatCh for MultiLocus Sequence Typing**
To run ReMatCh in a set of samples for MLST, either by providing the list of IDs/taxon name for download or the directory containing the sample forlders, the option `--mlst` needs to be provided with species name (same as in PubMLST) to be used in MLST determination. If more than one scheme is available for the species, the desired schema number should be passed to ReMatCh with the `--mlstSchemaNumber` option.
A fasta file containing the MLST reference sequences (`-r`) is required, with the size of the flanking regions  enough to allow the alignment of one read (set with the opion `--extraSeq`). In our experience, the addition of 200nt upstream and downstream of the target region when using Illumina Miseq data (150nt reads), will have the desired effect, and these flanking regions will be ignored in variant calling, unless there is an INDEL affecting the target sequence. Alternatively the `--mlstReference` option can be used, telling ReMatCH to use the curated scheme for the MLST scheme, if available, as reference sequences with 200nt flanking the target regions, or the first alleles of each MLST gene fragment in PubMLST as reference sequences.
As default, ReMatCh uses the consensus sequence "noMatter" in MLST determination, but this can be changed with the `--mlstConsensus` option. IF the option `--doubleRun` is used, ReMatCh can determine the MLST for the second run only, or for both runs, with the `--mlstRun` option. By default the MLST will be determined in both runs.
The MLST results will be in the *mlst_report.*.tab* in the `--workdir`.



Outputs
-------
**run.*.log**  
ReMatCh running log file.  

**sample_report.*.tab**
 - *sample* - Sample ID
 - *sample_run_successfully* - Reports whether the sample globally run successfully
 - *sample_run_time* - Global sample running time (in seconds)
 - *files_size* - Sum of files size (in bytes)
 - *download_run_successfully* - Reports whether the sample downloading (if requested) run successfully
 - *download_run_time* - Download running time
 - *rematch_run_successfully_first* - Reports whether the first run of ReMatCh module run successfully
 - *rematch_run_successfully_second* - Reports whether the second run of ReMatCh module run successfully
 - *rematch_run_time_first* - ReMatCh first running time
 - *rematch_run_time_second* - ReMatCh second running time
 - *number_absent_genes_first* - Number of absent genes determined in the first ReMatch run
 - *number_genes_multiple_alleles_first* - Number of genes with multiple alleles among the genes present determined in the first ReMatch run
 - *mean_sample_coverage_first* - Mean sample coverage depth (only considering the genes present) determined in the first ReMatch run
 - *number_absent_genes_second* - Number of absent genes determined in the second ReMatch run
 - *number_genes_multiple_alleles_second* - Number of genes with multiple alleles among the genes present determined in the second ReMatch run
 - *mean_sample_coverage_second* - Mean sample coverage depth (only considering the genes present) determined in the second ReMatch run
 - *run_accession* - ENA Run Accession number used to download
 - *instrument_platform* - Sequencing technology used reported by ENA
 - *instrument_model* - Instrument model used for sequencing reported by ENA
 - *library_layout* - Single or paired-end sequencing reported by ENA
 - *library_source* - Library type reported by ENA (e.g. genomic, transcriptomic, synthetic)
 - *extra_run_accession* - Extra ENA Run Accession numbers found for the ID provided for download
 - *date_download* - Date of downloading try
 - *fastq_used* - Fastq files used in ReMatCh module

**combined_report.data_by_gene.*.tab**  
_combined_report.data_by_gene.first_run.\*.tab_ and _combined_report.data_by_gene.second_run.\*.tab_  
This file contains a report with gene (in columns) presence/absence and coverage depth for the different samples (in lines).  
In the case of genes being present (genes with at least `--minGeneCoverage` percentage of target reference gene sequence covered with `--minCovPresence` reads and with at least `--minGeneIdentity` percentage identity of target reference gene sequence covered with `--minCovCall` reads), the script will provide the mean target sequence coverage, otherwise will report "absent_" for genes not present.  
In case of multiple alleles occurrence, if the frequency of the dominant allele is lower than `--minFrequencyDominantAllele` and the frequency of the most frequent minority allele is higher than 50% of the total of the minority alleles or is 50% but only 2 minority alleles exist, "multiAlleles_" will be reported.  

**cpu_information.*.cpu.txt** and **cpu_information.*.slurm.txt**  
Store CPUs and SLURM information at the time of run.  

**mlst_report.*.tab**
This file contains a report with the MLST information (columns) for the different samples (in lines).
For each sample, the file will have information on the run the MLST was determined (first or second), the consensus sequenced used (noMatter, correct or alignment), the ST obtained ( or '-' if no ST was obtained) and the allele number (or '-' if not an exact match) for each loci in the scheme.

**Samples folders**  
For each sample, three fasta files will be produced:  
 - *sample.noMatter.fasta* - Fasta file containing the target gene sequence with the more probable nucleotides (determined by the majority rule for positions covered by >= *--minCovPresence*). Positions with less than `--minCovPresence` coverage depth will be considered as deletions.
 - *sample.correct.fasta* - Fasta file containing the target gene sequence with the correct nucleotides. Positions with less than `--minCovPresence` coverage depth will be considered as deletions, with less than `--minCovCall` coverage depth will be coded as "N" (due to low certainty in calling SNP) and positions with possible multiple alleles will also be considered as "N".
 - *sample.alignment.fasta* - Will be exactly the same as *sample.correct.fasta*, but INDELs will be coded as "N" in order to produce sequences with the same size that can be then concatenated to produce an alignment file
 - *rematchModule_report.txt* - Report file containing gene information: 1) gene name, 2) percentage of target gene sequence covered with at least `--minCovPresence` read depth, 3) Mean target gene coverage depth of present positions, 4) percentage of target gene sequence with lower `--minCovCall` coverage depth, 5) number of positions in target gene sequence containing multiple alleles, 6) percentage identity of target gene sequence covered with at least `--minCovCall` read depth. The general sample information will also be stored: number of absent genes, number of genes with multiple alleles among the genes present and the mean sample coverage depth (only considering the genes present).
 - *rematch_module/* - Folder containing the temporary files. Only kept if `--debug` option is specified. It will contain the *alignment.bam*, bam and fasta indexes, *sequence_data/* folder with subfolders (named with numbers) for each sequence in `--reference` file. In each sequence folder the different consensus _\*.vcf_ files and the original _samtools_mpileup.\*.vcf_ and _samtools_depth.\*.vcf_ files
 - *rematch_second_run/* - Folder containing the same files/folders described above, but for the *second_run*. Only created if `--doubleRun` is set



Contact
-------
Miguel Machado  
<mpmachado@medicina.ulisboa.pt>

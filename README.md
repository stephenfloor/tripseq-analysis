# geseq-analysis
Analysis tools for GE-seq data.  There are various individual tools inside to perform different tasks.  Limited documentation below, contact floor@berkeley.edu with questions. 


## transcriptome_properties.py 

Calculate the properties of an input transcriptome (or regions thereof). Input format is BED, output files are .csv files with various properties as specified on the command line.

### usage

usage: transcriptome_properties.py [-h] -i INPUT -g GENOME [--gc] [--length]
                                   [--exonct] [--nt NT] [-o OUTPUT]
                                   [--window WINDOW]
                                   [--convtorefseq CONVTOREFSEQ]
                                   [--targetscanfile TARGETSCANFILE]
                                   [--deltag] [--lfold] [--cap-structure]
                                   [--kozak] [--uorf-count] [--uorf-overlap]
                                   [--start-codon] [--rare-codons]
                                   [--mirna-sites] [--au-elements]


optional arguments:
  -h, --help            show this help message and exit

Global arguments:
  -i INPUT, --input INPUT
                        The transcriptome region (BED format)
  -g GENOME, --genome GENOME
                        The genome for the input transcriptome
  --gc                  Calculate GC content
  --length              Calculate length
  --exonct              Count # of exons
  --nt NT               Number of threads (default is 8 or 4 for lfold)
  -o OUTPUT, --output OUTPUT
                        Output basename (e.g. CDS)
  --window WINDOW       Window size for sliding window calculations (default
                        75)
  --convtorefseq CONVTOREFSEQ
                        Filename to convert input annotations to refseq (for
                        targetscan; e.g. knownToRefSeq.txt)
  --targetscanfile TARGETSCANFILE
                        Filename of targetscan scores (e.g.
                        Summary_Counts.txt)
  --deltag              Calculate min deltaG in sliding window of size
                        --window over region
  --lfold               Use RNALfold to calculate MFE rather than RNAfold
                        (faster but does not compute centroid,MEA)

5' UTR specific arguments:
  --cap-structure       Calculate structure at the 5' end

Start-codon-specific arguments:
  --kozak               Calculate Kozak context score
  --uorf-count          Calculate number of 5' UTR uORFs (starting with
                        [ACT]TG)
  --uorf-overlap        Overlap of uORF with start codon (implies --uorf-
                        count)
  --start-codon         Record the start codon used (ATG or other)

CDS-specific arguments:
  --rare-codons         Calculate codon usage properties

3' UTR specific arguments:
  --mirna-sites         Compile miRNA binding site info from targetscan
  --au-elements         Count number of AU-rich elements in the 3' UTR


###Requirements: 
* ViennaRNA RNAfold and RNALfold (http://www.tbi.univie.ac.at/RNA)
* HumanCodonTable (this page)
* AnnotationConverter (this page)
* TargetscanScores (this page) 
* SNFUtils (this page) 

## compare_geseq_clusters.py 

Take two lists of GEseq data (i.e. clusters) and compare them for genes that have the transcripts in each of the two different sets.  For each set of gene-linked transcript isoforms, compare input transcriptome features as calculated using transcriptome_properties.py 

### Usage

usage: compare_geseq_clusters.py [-h] --set1 FNAME ID ... [FNAME ID ... ...]
                                 --set2 FNAME ID ... [FNAME ID ... ...]
                                 --tx-to-gene TX_TO_GENE [-o OUTPUT] -n NREP
                                 [--txome-props TXOME_PROPS [TXOME_PROPS ...]]
                                 [--control] --txome-gtf TXOME_GTF

optional arguments:
  -h, --help            show this help message and exit
  --set1 FNAME ID ... [FNAME ID ... ...]
                        Files and IDs containing transcript distributions;
                        compare between set1 and set2
  --set2 FNAME ID ... [FNAME ID ... ...]
                        Files and IDs containing transcript distributions;
                        compare between set1 and set2
  --tx-to-gene TX_TO_GENE
                        Mapping between transcript ID in input file and gene
                        ID
  -o OUTPUT, --output OUTPUT
                        Output filename (default is stdout)
  -n NREP, --nrep NREP  Number of replicates of each point
  --txome-props TXOME_PROPS [TXOME_PROPS ...]
                        List of files with transcriptome properties to
                        correlate among (wildcards ok)
  --control             Perform randomized comparisons of input transcripts as
                        a control.
  --txome-gtf TXOME_GTF
                        Path to transcriptome GTF

### Requirements

* GTF.py (this page)
* Transcript.py (this page)
* SNFUtils.py (this page) 
* Two lists of transcripts to compare (i.e. clusters) 
* Lists of transcriptome properties to compare between transcript isoforms of the same gene in the two sets (generated by transcriptome_properties.py) 
* File containing transcript ID to gene mapping

## plot_geseq_transcript.py 

Plot an individual transcript or all transcripts of a gene.  Requires input polysome sequencing data (i.e. GEseq) or some other distribution. 

### Usage

usage: plot_geseq_transcript.py [-h] -i INPUT [-o OUTPUT] -n NREP --id ID
                                --tx-to-gene TX_TO_GENE [--text]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        File containing transcript distributions
  -o OUTPUT, --output OUTPUT
                        Output filename (default is stdout)
  -n NREP, --nrep NREP  Number of replicates of each point
  --id ID               Transcript ID(s) to print (can be partial; can be
                        comma-separated list)
  --tx-to-gene TX_TO_GENE
                        File containing transcript ID to gene name mapping
  --text                Output text data in addition to plots.

### Requirements

* SNFUtils.py (this page) 
* Input per-transcript distributions
* File containing transcript ID to gene mapping (if per-gene plotting is desired) 
  
## fpkm_to_tpm.py

Converts between FPKM and TPM (transcripts per million).  Uses the formula TPM_i = FPKM_i * 1e6 / sum(FPKM_g for all genes g)
Citation: http://lynchlab.uchicago.edu/publications/Wagner,%20Kin,%20and%20Lynch%20%282012%29.pdf

### Usage

usage: fpkm_to_tpm.py [-h] -i INPUT [-t SEPARATOR] [-o [OUTPUT]]
                      [--ignore IGNORE] [--filter FILTER] [-u]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        File containing identifiers to use for the merge
  -t SEPARATOR, --separator SEPARATOR
                        Field separator (default comma; "tab" for tabs;
                        "space" for whitespace
  -o [OUTPUT], --output [OUTPUT]
                        File to output to (default stdout)
  --ignore IGNORE       Number of columns to ignore (one-based; 1 ignores the
                        first column)
  --filter FILTER       Filter genes with TPM below arg
  -u, --unique          Only output lines with unique entries in column 1

### Requirements 
* A file with FPKM values to convert to TPM 

## Utility classes

### AnnotationConverter.py 

A class to provide for conversion between two annotation sets. 

### GTF.py 

A class to read GTF files - downloaded from https://gist.github.com/slowkow/8101481 and minimally modified 

### SNFUtils.py

A file providing various utility functions.

### HumanCodonTable.py

A class harboring information on human codon usage.

### TargetscanScores.py

A class to read in targetscan scores and provide accessor functions. 

### Transcript.py

A class defining a transcript and structural features associated with it. 

## compare_geseq_clusters.py 

### Usage

### Requirements

## AnnotationConverter.py  GTF.py  HumanCodonTable.py   SNFUtils.py  TargetscanScores.py  Transcript.py  fpkm_to_tpm.py 


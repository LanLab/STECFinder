# STECFinder
Clustering and Serotyping of Shigatoxin producing E. coli (STEC) using genomic cluster specific markers.

This tool can identify the serotype of STEC using cluster-specific genes and O-antigen/H-antigen genes. Input is either
illumina reads (fastq.gz) or genome assemblies (fasta).

---
# Dependencies
1. python (v3.6 or greater)
2. kma  (v1.3.15 or greater)
3. blastn (v2.9 or greater)
---

# Installation 
## Option 1: Clone repository from gitHub
````
git clone https://github.com/LanLab/STECFinder.git

cd STECFinder

python setup.py install
````
Make sure that you have the dependencies installed.

## Option 2: Conda installation
`conda install -c bioconda -c conda-forge stecfinder`

[![Anaconda-Server Badge](https://anaconda.org/bioconda/stecfinder/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/stecfinder/badges/downloads.svg)](https://anaconda.org/bioconda/stecfinder)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/stecfinder/badges/version.svg)](https://anaconda.org/bioconda/stecfinder)

# Usage:
```commandline
usage: 
STECFinder.py -i <input_data1> <input_data2> ... OR
STECFinder.py -i <directory/*> OR 
STECFinder.py -i <Read1> <Read2> -r [Raw Reads]

Input/Output:
  -i I [I ...]          <string>: path/to/input_data (default: None)
  -r                    Add flag if file is raw reads. (default: False)
  -t T                  number of threads. Default 4. (default: 4)
  --hits                shows detailed gene search results (default: False)
  --output OUTPUT       output file to write to (if not used writes to stdout and tmp folder in current dir) (default: None)

Misc:
  -h, --help            show this help message and exit
  --check               check dependencies are installed (default: False)
  -v, --version         Print version number (default: False)

Algorithm cutoffs:
  --cutoff CUTOFF       minimum read coverage for gene to be called (default: 10.0)
  --length LENGTH       percentage of gene length needed for positive call (default: 50.0)
  --ipaH_length IPAH_LENGTH
                        percentage of ipaH gene length needed for positive gene call (default: 10.0)
  --ipaH_depth IPAH_DEPTH
                        When using reads as input the minimum depth percentage relative to genome average for positive ipaH gene call (default: 1.0)
  --stx_length STX_LENGTH
                        percentage of stx gene length needed for positive gene call (default: 10.0)
  --stx_depth STX_DEPTH
                        When using reads as input the minimum depth percentage relative to genome average for positive stx gene call (default: 1.0)
  --o_length O_LENGTH   percentage of wz_ gene length needed for positive call (default: 60.0)
  --o_depth O_DEPTH     When using reads as input the minimum depth percentage relative to genome average for positive wz_ gene call (default: 1.0)
  --h_length H_LENGTH   percentage of fliC gene length needed for positive call (default: 60.0)
  --h_depth H_DEPTH     When using reads as input the minimum depth percentage relative to genome average for positive fliC gene call (default: 1.0)
```


# Example:

Run on a folder containing pairs of fastq files using kma for gene identification

```commandline
python STECfinder.py -r -i "/input/reads/folder/*" --output "/output/file/name"
```
Run on a folder containing genome files using kma for gene identification
```commandline
python STECfinder.py -i "/input/genomes/folder/*" --output "/output/file/name"
```


# Output:
````
Sample	        Cluster	Cluster Serotype	Serotype	Big10 serotype	O antigens	H antigens	stx type	ipaH presence	Notes
SRR3995879	O157H7	O157:H7	                O157:H7	        -	        wzy_O157	H7	        stx2a	        -	        -
SRR1917514	C4	O5:H9	                O5:H9	        -	        wzy_O5	        H9	        stx1a	        -	        -
````

## Column descriptions:
| Column | Description |
| ----------- | ----------- |
|Sample|Input strain ID (extracted from input files)|
|Cluster| Phylogenetic STEC cluster predicted by accessory genome specific gene sets|
|Cluster Serotype| Filtered serotype from antigen gene matches but only allowing serotypes previously seen in the predicted cluster|
|Serotype| unfiltered serotype from antigen gene matches|
|Big10 serotype| The serotype of the isolate as predicted by accessory genome serotype specific gene sets for the top 10 non O157:H7 STEC serotypes|
|O antigens| o antigen gene matches|
|H antigens| h antigen gene matches|
|_stx_ type| _stx_ toxin genes detected|
|_ipaH_| presence or absence of _ipaH_ gene|
|Notes| Other information if unexpected results are observed|

### Note on stx2 allele call outputs
- \* denotes a call with some uncertainty (either a non perfect match to a known allele for stx1, or a minority of allele specific SNPs for stx2)
- calls separated by "/" are multiple possible alleles for a single stx2 locus
- calls separated by "," are calls for multiple separate stx2 loci in the same genome

Assemblies will often merge multiple stx2 genes into one assembled locus. Therefore it is only possible to detect multiple stx2 alleles in one strain using raw reads, which allow the frequencies of stx2 type defining SNPs to be evaluated.

## Column descriptions for additional tables produced by --hits flag:
### - PROCESSED GENE SET -
Filtered gene hits used to make clustering, serotyping and gene presence calls for algorithm

| Column | Description |
| ----------- | ----------- |
|gene|Locus ID in resources/genes.fasta|
|subject_length| Length of locus|
|perc_ident| Percentage identity of matching region|
|length_percentage| percentage of the length of the subject sequence that is matched to|
|score| blast: bitscore, KMA: mapping score|
|gene_type| category of locus|
|depth|(-r, read input only) depth of reads matching to gene|
|normalised_depth|(-r, read input only) depth but as a percentage of average depth of 7 housekeeping genes|

### -RAW KMA HITS- and -RAW BLAST HITS- 
Raw output tables for gene matching programs. Columns as per documentation of those tools

BLAST:

| Column | Description |
| ----------- | ----------- |
|sseqid | Subject Seq-id|
|slen | Subject sequence length|
|length | Alignment length|
|sstart | Start of alignment in subject|
|send | End of alignment in subject|
|pident | Percentage of identical matches|
|bitscore | Bit score|

[KMA](https://bitbucket.org/genomicepidemiology/kma/src/master/KMAspecification.pdf):

| Column | Description |
| ----------- | ----------- |
|Template| Contains the name of the template, default is the fasta header from the template sequence, including any spaces, tabs or special characters.
|Score| Is the ConClave score (accumulated alignment score), from all reads that were accepted tomatch this template.|
|Expected| Is the expected Score, if all mapping reads were normally distributed over the entire database.|
|Template_length| Is the length of the template sequence, without preceding and trailing N’s.|
|Template_Identity| Is the number of bases in the consensus sequence that are identical to the template sequence divided by the Template_length. In other words, the percentage of identical nucleotides between template and consensus w.r.t. the template.|
|Template_Coverage| Is the percentage of bases in the template that is covered by the consensus sequence. A Template_Coverage above 100% indicates the presence of more insertions than deletions.|
|Query_Identity| Is the number of bases in the template sequence that are identical to the consensus sequence divided by the length of the consensus. In other words, the percentage ofidentical nucleotides between template and consensus w.r.t. the consensus.|
|Query_Coverage| Are the reciprocal values of the Template_Coverage. A Query_Coverage above 100% indicates the presence of more deletions than insertions.
|Depth| Is the depth of coverage of the template. Commonly referred to as X-coverage, coverage, abundance, etc.|
|Q_value| Is the obtained quantile in a --- distribution, when comparing the obtained Score with the Expected, using a McNemar test.|
|P_value| Is the obtained p-value from the quantile Q_value.|

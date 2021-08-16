# STECFinder
Clustering and Serotyping of Shigatoxin producing E. coli using genomic cluster specific markers
This is a tool that is used to identify the serotype of STEC using cluster-specific genes and O-antigen/H-antigen genes.

---
# Dependencies
1. samtools (v1.10)
2. python (v3.7.3)
3. bwa (v0.7.17-r1188)
4. kma  (v1.3.15)
---
# Installation 

install shigeifinder from conda: 
`conda install -c bioconda -c conda-forge shigeifinder`

install kma in conda: 
`conda install -c bioconda kma`

Usage:
```commandline
 
STECFinder.py -i <input_data1> <input_data2> ... OR
STECFinder.py -i <directory/*> OR 
STECFinder.py -i <Read1> <Read2> -r [Raw Reads]

optional arguments:
  -h, --help       show this help message and exit
  -i I [I ...]     <string>: path/to/input_data
  -r               Add flag if file is raw reads.
  -t T             number of threads. Default 4.
  --hits           To show the blast/alignment hits
  --dratio         To show the depth ratios of cluster-specific genes to House Keeping genes
  --update_db      Add flag if you added new sequences to genes database.
  --output OUTPUT  output file to write to (if not used writes to stdout and tmp folder in current dir)
  --check          check dependencies are installed
```

#Example:

Run on a folder containing pairs of fastq files using kma for gene identification

```commandline
python STECfinder.py -r -i "/input/reads/folder/*" --output "/output/file/name"
```
Run on a folder containing genome files using kma for gene identification
```commandline
python STECfinder.py -i "/input/genomes/folder/*" --output "/output/file/name"
```

#Output:
````
#SAMPLE	        STX	cluster	big10_serotype	serotype	cluster_serotype    oantigens	hantigens	IPAH	Notes
SRR3995879	stx2a	O157H7  -               O157:H7	        O157:H7	            wzy_O157	H7	        -	-
SRR1917514	stx1a	C4      -               O5:H9	        O5:H9	            wzy_O5	H9	        -	-
````
###Column descriptions:
SAMPLE: Input strain ID 
STX:
cluster
big10_serotype
serotype
cluster_serotype
oantigens
hantigens
IPAH
Notes
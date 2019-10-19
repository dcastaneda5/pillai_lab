# Amplicon Deep Sequencing Compared to Traditional msp genotyping of *Plasmodium falciparum* malaria

## The contents of this folder are related to the data used for the research article "Amplicon deep sequencing compared to traditional msp genotyping of P. falciparum malaria" by Pillai lab.

### The software and commands used for this project include:

SRST2 (to confirm the presence of alleles in the built database):
```
srst2 --input_pe <sample1_f.fastq> <sample1_r.fastq> --min_coverage 90 --gene_db <srst2_whole_db.fasta> --threads 16
```

VSEARCH (to linearize, dereplicate, and add abundance values):
```
vsearch --fastq_filter <sample1.fastq> --fastaout <sample1.fasta>

vsearch --derep_fulllength <sample1.fasta> --sizeout --relabel_sha1 --fasta_width 0 --output <sample1_linearized.fasta>
```

SWARM V.2 (to create and visualize the clustering output):
```
swarm -d 2 -t 8 -z -i <internal_file> -w <otu_file> <sample1_linearized.fasta> > <sample1_out.swarms>

python graph_plot.py -i <internal_file> -s <sample1_out.swarms> -d 0
```

RStudio (to perform the wilcoxon test and visualize the significance of the results):

Three R files were used for statistical analysis, which are found in the scripts folder.

Basic stats:

Basic descriptive statistics was achieved by using the `basic_stats_fastq.py` in the scripts folder.





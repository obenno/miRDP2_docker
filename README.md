# Docker container recipe of miRDP2 pipeline for miRNA identification

## miRDP2

The miRPD2 is a flexible pipeline for identification of miRNA from
miRNA-seq data, with updated criteria for plant. The source code
was hosted on [sourceforge](https://sourceforge.net/projects/mirdp2/).
This images was built based on miRDP2 v1.1.4, check the detail from
it's [manual](https://sourceforge.net/projects/mirdp2/files/version%201.1.4/miRDP2_manual-v1.1.4.pdf/download).

### Cite

Zheng Kuang, Ying Wang, Lei Li, Xiaozeng Yang, miRDeep-P2: 
accurate and fast analysis of the microRNA transcriptome in 
plants, Bioinformatics, Volume 35, Issue 14, July 2019, 
Pages 2521–2522, https://doi.org/10.1093/bioinformatics/bty972

## Data preparation

### ncRNA contaminant

The input reads will have to be filtered to discard contaminants
from other short RNAs (tRNA, snRNA, snoRNA, ribozyme etc.) before
being mapped to the genome.


The ncRNAs other than miRNA were extracted from [Rfam](https://rfam.xfam.org/)'s datasets.
The family information and fasta files were downloaded from Rfam's FTP.
server.

```
## Download ncRNA fasta
wget -c http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz
zcat family.txt.gz |cut -f 1,19 | awk -F"\t" '$2~/tRNA|snRNA|snoRNA|ribozyme/'| awk '{print "wget -c --directory-prefix=./fasta/ http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/"$1".fa.gz"}' | bash

## Build bowtie index
cat fasta/RF*.fa.gz | gzip -c > miRDP2_ncRNA_rfam_v14.5.fa
bowtie-build -f ./miRDP2_ncRNA_rfam_v14.5.fa /home/ubuntu/Tools/miniconda3/envs/mirdeep-p2/bin/scripts/index/rfam_index
```





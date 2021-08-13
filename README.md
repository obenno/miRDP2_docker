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

## Image

docker pull obenno/mirdp2_pipeline

## Data preparation

### ncRNA contaminant

The input reads will have to be filtered to discard contaminants
from other short RNAs (tRNA, snRNA, snoRNA, ribozyme etc.) before
being mapped to the genome.


The ncRNAs other than miRNA were extracted from [Rfam](https://rfam.xfam.org/)'s datasets (v14.5).
The family information and fasta files were downloaded from Rfam's FTP server.

```
## Download ncRNA fasta
wget -c http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz
zcat family.txt.gz |cut -f 1,19 | awk -F"\t" '$2~/tRNA|snRNA|snoRNA|ribozyme/'| awk '{print "wget -c --directory-prefix=./fasta/ http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/"$1".fa.gz"}' | bash

## Build bowtie index
cat fasta/RF*.fa.gz | gunzip -c > miRDP2_ncRNA_rfam_v14.5.fa
bowtie-build -f ./miRDP2_ncRNA_rfam_v14.5.fa ./index/rfam_index
```

The `rfam_index` was prepared and packaged into the image.

### Known miRNA 

If refFile provided, miRNA will be compared with reference set, and be named accordingly
The refFile should contain 6 columns: 5p_seq, 3p_seq, 5p_name, 3p_name, precursor_name, precursor_seq
all known miRNA should be sorted and uniqufied by 5p_seq, 3p_seq (bedtools groupby)

example for prepare refFile from [PmiREN2.0 soybean known miRs](https://www.pmiren.com/ftp-download/Glycine_max_Gma/Glycine_max_basicInfo.txt)

```
cut -f 1,13,15,18,19,22 Glycine_max_basicInfo.txt |
    awk '
    BEGIN{OFS="\t"}
    NR==1{print}
    NR>1{k=length($2); sub("*","", $5); if(index($2,$4)<=k/2){$3=$3"-5p";$5=$5"-3p"}else{$3=$3"-3p";$5=$5"-5p"}; print}
    ' |
    awk 'NR>1{if($3~/5p$/){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}else{print $1"\t"$2"\t"$5"\t"$6"\t"$3"\t"$4}}' |
    awk '{print $3"\t"$4"\t"$5"\t"$6"\t"$1"\t"$2}' |
    sort -k 2,2 -k 4,4 -k 1,1 -k 3,3 |
    bedtools groupby -g 2,4 -c 1,3,5,6 -o collapse > mature_precursor_seq.tsv
```

## Running Command

Assume the input files are in the current directory, the reference genome and its bowtie idnex file are
located in `/home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/assembly`. The current directory will be 
bound to `/data`, and folder containing reference genome and index will be bound to `/database`. Reference
genome file is `Gmax_508_v4.0.fa`, and its index name is also `Gmax_508_v4.0.fa`.

```
docker run -it -d --rm -v /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/assembly:/database \
           -v $(pwd):/data -u $(echo $UID):1000 obenno/mirdp2_pipeline:0.01 \
           -g /database/Gmax_508_v4.0.fa \
           -x /database/Gmax_508_v4.0.fa \
           -b /data/test_inputList \
           -t -d /data/mature_precursor_seq.tsv -p 40 -o /data/OutputDir
```

## Full Options

Main script is modified from miRDP2 v1.1.4.

```
miRDP2-v1.1.4_pipeline.bash <OPTIONS>

PLEASE REFER TO MIRDP2 MANUAL ALSO.

If input is single file with fasta format, then it will be treated
as pre-formatted fasta, and no trimming or collapsing will be conducted.

  OPTIONS:
    NECESSARY:
    -g/--genome <file>        reference genome file in fasta format
    -x/--index <path>/prefix  bowtie index file corresponding to reference genome
    -f/--fasta                formatted sRNA-seq file in fasta format
    -q/--fastq                unformatted sRNA-seq file in fastq
    -i/--input                input file, use batch option for multiple input files
    -b/--batch                list file of input files
    -o/--output <path>        path to output folder

    OPTIONAL:
    -d/--db                   reference miRNA information in tsv format (for family identification)
    -t/--trim                 switch for adapter trimmer step, default false
    -a/--adapter              adapter sequence to be trimmed [TGGAATTCTCGGGTGCCAAGG]
    -L/--locate <int>         THRESHOLD of different locations one read can map to, default is 15
    -N/--length <int>         THRESHOLD of length of precursor candidates, default is 300
    -M/--mismatch <int>       THRESHOLD of allowed mismatches for bowtie mapping, default is 0
    -R/--rpm <int>            THRESHOLD of reads rpm when filtering, default is 10
    -p/--thread <int>         number of thread used for RNAfold, default is 1
    --large-index             use this option when using large bowtie index, e.g. when using *.ebtwl files
    -h/--help                 print this help

```

## Dockerfile

Dockerfile was written referring to https://pythonspeed.com/articles/activate-conda-dockerfile/
and Dockstore's [docs](https://docs.dockstore.org/en/stable/getting-started/getting-started-with-docker.html).

## Issues

The `-j` of `RNAfold` doesn't work in the container.


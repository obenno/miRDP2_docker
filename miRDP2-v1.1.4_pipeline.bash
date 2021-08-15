#!/bin/bash

# The main pipeline of miRDP2


usage="$0 <OPTIONS>

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
"

	
#####################  default parameters  #####################
results_folder="."
locate=15
length=300
mismatch=0
rpm=10
large=""
thread=1
bt2tag="false"

tag="f"
batch=""
trim="false"
refFile=""
adapter="TGGAATTCTCGGGTGCCAAGG"

#####################  get input  #####################

if [[ -z $1 ]]; then
    echo -e "$usage";
    exit 1;
fi;

while [[ ! -z $1 ]]
do
  case "$1" in
	-g|--genome)  genome=$2; shift 2;;	#input genome file
    -x|--index)   bowtie_index=$2; shift 2;;	#input index file
    -f|--fasta)   tag="f"; shift;;	#input fasta seq file
	-q|--fastq)   tag="q"; shift;;  #input fastq seq file
    -t|--trim)    trim="true"; shift;;
    -a|--adapter) adapter=$2; shift;;
	-i|--input)   input=$2; shift 2;;  #input file(s)
	-b|--batch)   batch=$2; shift 2;;  #input file list file
    -o|--output)  results_folder=$2; shift 2;;	#output folder; default is `.'
    -L|--locate)  locate=$2; shift 2;;	#multi_location threshold
    -N|--length)   length=$2; shift 2;;	#pre length
	-M|--mismatch) mismatch=$2; shift 2;;	#mismatch
    -R|--rpm)      rpm=$2; shift 2;;	#RPM threshold
    -p|--thread)   thread=$2; shift 2;;	#thread number
	##-T|--bowtie2)	  bt2tag="BT2"; shift;;
    -d|--db)       refFile=$2; shift 2;;
	--large-index)    large="--large-index"; shift;;
    -h|--help)    echo -e "$usage"; exit 1;;
    *)            echo -e "$usage"; exit 1;;
  esac
done



#####################  change variables  #####################

##src=${0%%/miRDP2-v1.1.4_pipeline.bash}
## Set src path
src="/opt/conda/envs/miRDP2/bin"
##src="/home/ubuntu/Tools/miniconda3/envs/miRDP2/bin"

mirdp=$src/scripts			#miRDP script folder
len=$locate				#allowed maximum mapping site number
var=$mismatch				#allowed mismatch number
threshold=$rpm		#threshold for filter reads (rpm)

#filename=$filename-$len-$var-$threshold

## Print arguments
## echo "\$input: $input"
## echo "\$batch: $batch"
## echo "\$tag: $tag"


#####################  begin pipeline  #####################
## Original pipeline will process each input file independently
## This modified version will merage all the input fastq into one
## before identifying miRNAs
filename="merged_out"
mkdir -p $results_folder/$filename
echo "begin pipeline" $(date +[\ %F\ %X\ ]) > $results_folder/${filename}/progress_log
echo "" > $results_folder/${filename}/script_err
echo "" > $results_folder/${filename}/script_log

# preprocess of fastq file
## if [ $tag = "q" ]; then
## 	perl -e '$c=1; while(<>){if(/^@/){$seq=<>;chomp($seq); <>;<>; $hash{$seq}++}} foreach my $s(sort keys %hash){printf ">read%08d_x$hash{$s}\n$s\n",$c; $c++;}' $seq > $results_folder/${filename}/${filename}.ori.fa
## 	seq="$results_folder/${filename}.ori.fa"
## fi

# preprocess & mapping
## bowtie2 disabled
## if [ $bt2tag == "BT2" ]; then
## 	#filter reads -- mapping to ncRNA seq (rRNA, tRNA, etc)
## 	bowtie2 -x $mirdp/index/rfam_index -f $seq > $results_folder/${filename}/rfam_reads.sam 2>> $results_folder/${filename}/script_err
## 	#filter reads -- mapping to known miR mature seq
## 	bowtie2 -x $mirdp/index/mature_index -f $seq > $results_folder/${filename}/known_miR.sam 2>> $results_folder/${filename}/script_err
## 	#filter reads -- extract reads with high count number
## 	perl $mirdp/preprocess_reads-SAM.pl $seq $results_folder/${filename}/rfam_reads.sam $results_folder/${filename}/known_miR.sam $threshold $results_folder/${filename}/${filename}.fa $results_folder/${filename}/${filename}-processed.fa $results_folder/${filename}/${filename}.total_reads 2>> $results_folder/${filename}/${filename}_err
## 	
## 	#mapping filtered reads
## 	bowtie2 -a -N 1 $large -x $bowtie_index -f $results_folder/${filename}/${filename}-processed.fa > $results_folder/${filename}/${filename}_processed.sam 2>> $results_folder/${filename}/script_err
## 	#filter reads -- ignore reads mapping to too many sites
## 	perl $mirdp/convert_SAM_to_blast.pl $results_folder/${filename}/${filename}_processed.sam $results_folder/${filename}/${filename}.fa $genome $var > $results_folder/${filename}/${filename}-processed.bst 2>>$results_folder/${filename}/${filename}_err
## fi

## Test input and batch option, at least one of the two is required
if [[ ! -z $input ]] || [[ -f $batch ]]; then
    if [[ $input =~ .(fa|fasta).gz$ && tag == "f" ]]; then
        formattedInput=$(mktemp -p $results_folder "formatted.XXXXXXXXXX")
        zcat $input > $formattedInput
    elif [[ $input =~ .(fa|fasta)$ && $tag == "f" ]]; then
        formattedInput=$(mktemp -p $results_folder "formatted.XXXXXXXXXX")
        cat $input > $formattedInput
    elif [[ $input =~ .(fq|fastq).gz$ && $tag == "q" ]]; then
        formattedInput=$(mktemp -p $results_folder "formatted.XXXXXXXXXX")
        if [[ $trim == "true" ]]; then
            inputTrimmed=$(mktemp -p $results_folder ${input%%.[fq|fa]*}".trimmed.XXXXXXXXXX")
            cutadapt -a $adapter -j $thread --fasta -o $inputTrimmed $input >> $results_folder/${filename}/script_err
            fastx_collapser -Q33 -i $inputTrimmed |
                awk '{if($1~/^>/){id=substr($1,2); split(id, tmp, "-"); print ">read"tmp[1]"_x"tmp[2]}else{print}}' > $formattedInput
        else
            zcat $input | fastx_collapser -Q33 |
                awk '{if($1~/^>/){id=substr($1,2); split(id, tmp, "-"); print ">read"tmp[1]"_x"tmp[2]}else{print}}' > $formattedInput
        fi
    elif [[ $input =~ .(fq|fastq)$ && $tag == "q" ]]; then
        formattedInput=$(mktemp -p $results_folder "formatted.XXXXXXXXXX")
        if [[ $trim == "true" ]]; then
            inputTrimmed=$(mktemp -p $results_folder ${input%%.[fq|fa]*}".trimmed.XXXXXXXXXX")
            cutadapt -a $adapter -j $thread --fasta -o $inputTrimmed $input >> $results_folder/${filename}/script_err
            fastx_collapser -Q33 -i $inputTrimmed |
                awk '{if($1~/^>/){id=substr($1,2); split(id, tmp, "-"); print ">read"tmp[1]"_x"tmp[2]}else{print}}' > $formattedInput
        else
            cat $input | fastx_collapser -Q33 |
                awk '{if($1~/^>/){id=substr($1,2); split(id, tmp, "-"); print ">read"tmp[1]"_x"tmp[2]}else{print}}' > $formattedInput
        fi
    elif [[ -f $batch ]]; then
        formattedInput=$(mktemp -p $results_folder "formatted.XXXXXXXXXX")
        if [[ $trim == "true" ]]; then
            batchTrimmed=$(mktemp -p $results_folder "batchTrimmed.XXXXXXXXXX")
            for i in $(cat $batch);
            do
                tmpTrimmed=$(mktemp -p $results_folder ${i%%.[fq|fa]*}".trimmed.XXXXXXXXXX")
                cutadapt -a $adapter -j $thread --fasta -o $tmpTrimmed $i >> $results_folder/${filename}/script_err
                echo $tmpTrimmed >> $batchTrimmed
            done
            awk 'BEGIN{a=""}{a=a" "$1}END{print "cat "a}' $batchTrimmed |
                bash |
                fastx_collapser -Q33 |
                awk '{if($1~/^>/){id=substr($1,2); split(id, tmp, "-"); print ">read"tmp[1]"_x"tmp[2]}else{print}}' > $formattedInput
        else
            awk 'BEGIN{gz=0;a=""}{if($1~/.gz$/){gz=1}; a=a" "$1}END{if(gz==0){print "cat "a}else{print "zcat "a}}' $batch |
                bash |
	            fastx_collapser -Q33 |
                awk '{if($1~/^>/){id=substr($1,2); split(id, tmp, "-"); print ">read"tmp[1]"_x"tmp[2]}else{print}}' > $formattedInput
        fi
    else
        echo "Either '--input' or '--batch' option is not valid."
        exit 1
    fi
else
    echo "Please provide input file with '--input' or '--batch' option."
    exit 1
fi


## Prefilter input reads to exclude ones with too few abundance
formattedInput_filtered=$(mktemp -p $results_folder "formattedFiltered.XXXXXXXXXX")
fasta_formatter -t -i $formattedInput |
    awk '{quant=substr($1, index($1, "_x")+2)+0; if(quant>=5){print ">"$1"\n"$2}}' > $formattedInput_filtered

if [ $bt2tag == "false" ]; then
    #filter reads -- mapping to ncRNA seq (rRNA, tRNA, etc)
    echo "Mapping reads to rfam ncRNAs (rRNA, tRNA, snRNA, snoRNA): " >> $results_folder/${filename}/script_err
    echo "bowtie -v 0 -x $mirdp/index/rfam_index -p $thread -f $formattedInput_filtered" >> $results_folder/${filename}/script_err
    bowtie -v 0 \
           -x $mirdp/index/rfam_index \
           -p $thread \
           -f $formattedInput_filtered \
           > $results_folder/${filename}/rfam_reads.aln \
           2>> $results_folder/${filename}/script_err
    #filter reads -- mapping to known miR mature seq
    echo "Mapping reads to miRBase: " >> $results_folder/${filename}/script_err
    echo "bowtie -v 1 -x $mirdp/index/mature_index -p $thread -f $formattedInput_filtered" >> $results_folder/${filename}/script_err
    bowtie -v 1 \
           -x $mirdp/index/mature_index \
           -p $thread \
           -f $formattedInput_filtered > $results_folder/${filename}/known_miR.aln 2>> $results_folder/${filename}/script_err
    #filter reads -- extract reads with high count number

    ## Discard reads with hit(s) of rfam_ncRNA dataset, and the count<10
    ## awk '{print}' $results_folder/${filename}/rfam_reads.aln |
    ##     sort -u |
    ##     awk '
    ##     NR==FNR{a[$1]}
    ##     NR>FNR{
    ##         if(!($1 in a)&&length($2)>=19&&length($2)<=24){
    ##             split($1, tmp, "_x");
    ##             count=tmp[2];
    ##             if(count>=10){
    ##                 print ">"$1"\n"$2
    ##             }
    ##         }
    ##     }' - <(fasta_formatter -t -i $seq) > $results_folder/${filename}/${filename}-processed.fa
    ## 
    ## fasta_formatter -t -i $results_folder/{filename}/${filename}-processed.fa |
    ##     awk 'BEGIN{a=0}{split($1,tmp,"_x"); a+=tmp[2]}END{print a}' > $results_folder/${filename}/${filename}.total_reads

    perl $mirdp/preprocess_reads.pl $formattedInput_filtered $results_folder/${filename}/rfam_reads.aln \
         $results_folder/${filename}/known_miR.aln \
         $threshold \
         $results_folder/${filename}/${filename}.fa \
         $results_folder/${filename}/${filename}-processed.fa \
         $results_folder/${filename}/${filename}.total_reads \
         2>> $results_folder/${filename}/${filename}_err
    ## $results_folder/${filename}/${filename}-processed.fa is filtered by
    ## 1. ncRNA
    ## 2. length (18, 24)
    ## 3. mapped to known mature.fa from mirbase or CPM > criteria
    ## will be used for excise precursor

    ## $results_folder/${filename}/${filename}.fa is only filtered by the
    ## first two criterias, will be used as signature and be mapped to precursor

    #mapping filtered reads
    echo "Mapping filtered reads to genome to excise precursors:" >> $results_folder/${filename}/script_err
    echo "bowtie -a -v $var $large -x $bowtie_index -p $thread -f $results_folder/${filename}/${filename}-processed.fa" >> $results_folder/${filename}/script_err
    bowtie -a -v $var $large \
           -x $bowtie_index \
           -p $thread \
           -f $results_folder/${filename}/${filename}-processed.fa \
           > $results_folder/${filename}/${filename}_processed.aln \
           2>> $results_folder/${filename}/script_err
    #filter reads -- ignore reads mapping to too many sites
    perl $mirdp/convert_bowtie_to_blast.pl \
         $results_folder/${filename}/${filename}_processed.aln \
         $results_folder/${filename}/${filename}.fa \
         $genome \
         > $results_folder/${filename}/${filename}-processed.bst \
         2>>$results_folder/${filename}/${filename}_err
fi

perl $mirdp/filter_alignments.pl \
     $results_folder/${filename}/${filename}-processed.bst \
     -c $len \
     > $results_folder/${filename}/${filename}-processed_filter${len}.bst \
     2>>$results_folder/${filename}/${filename}_err
echo "finish filtering" $(date +[\ %F\ %X\ ]) >> $results_folder/${filename}/progress_log


#excise precursor candidates and predict secondary strucure
perl $mirdp/excise_candidate.pl $genome $results_folder/${filename}/${filename}-processed_filter${len}.bst $length > $results_folder/${filename}/${filename}_precursors.fa 2>>$results_folder/${filename}/${filename}_err;

if [ "$thread" -gt 1 ]
then
	RNAfold --noPS -j$thread < $results_folder/${filename}/${filename}_precursors.fa > $results_folder/${filename}/${filename}_structures 2>>$results_folder/${filename}/${filename}_err;
else
	RNAfold --noPS < $results_folder/${filename}/${filename}_precursors.fa > $results_folder/${filename}/${filename}_structures 2>>$results_folder/${filename}/${filename}_err;
fi

echo "finish folding candidate precursors" $(date +[\ %F\ %X\ ]) >> $results_folder/${filename}/progress_log

#extract reads with no ncRNA for signature preparation
## if [ $bt2tag == "BT2" ]; then
## 	bowtie2 -a -N $var -a $large -x $bowtie_index -f $results_folder/${filename}/${filename}.fa > $results_folder/${filename}/${filename}.sam 2>> $results_folder/${filename}/script_err
## 	perl $mirdp/convert_SAM_to_blast.pl $results_folder/${filename}/${filename}.sam $results_folder/${filename}/${filename}.fa $genome $var > $results_folder/${filename}/${filename}.bst 2>>$results_folder/${filename}/${filename}_err
## fi

if [ $bt2tag == "false" ]; then
    echo "Mapping reads (non rfam_ncRNA) to the genome as candidate signatures:" >> $results_folder/${filename}/script_err
    echo "bowtie -a -v $var $large -x $bowtie_index -f $results_folder/${filename}/${filename}.fa" >> $results_folder/${filename}/script_err
    bowtie -a -v $var $large \
           -x $bowtie_index \
           -f $results_folder/${filename}/${filename}.fa \
           -p $thread \
           > $results_folder/${filename}/${filename}.aln \
           2>> $results_folder/${filename}/script_err
    perl $mirdp/convert_bowtie_to_blast.pl \
         $results_folder/${filename}/${filename}.aln \
         $results_folder/${filename}/${filename}.fa \
         $genome \
         > $results_folder/${filename}/${filename}.bst \
         2>>$results_folder/${filename}/${filename}_err
fi

perl $mirdp/filter_alignments.pl $results_folder/${filename}/${filename}.bst -c $len > $results_folder/${filename}/${filename}_filter${len}.bst 2>>$results_folder/${filename}/${filename}_err
perl $mirdp/filter_alignments.pl $results_folder/${filename}/${filename}_filter${len}.bst -b $results_folder/${filename}/${filename}.fa > $results_folder/${filename}/${filename}_filtered.fa 2>>$results_folder/${filename}/${filename}_err


#prepare reads signature file
if [[ ! -d $results_folder/${filename}/index ]]; then
    mkdir $results_folder/${filename}/index
fi
## if [ $bt2tag == "BT2" ]; then
## 	bowtie2-build -f $results_folder/${filename}/${filename}_precursors.fa $results_folder/${filename}/index/${filename}_precursors >> $results_folder/${filename}/script_log 2>> $results_folder/${filename}/script_err
## 	bowtie2 -a -x $results_folder/${filename}/index/${filename}_precursors -f $results_folder/${filename}/${filename}_filtered.fa > $results_folder/${filename}/${filename}_precursors.sam 2>> $results_folder/${filename}/script_err
## 	perl $mirdp/convert_SAM_to_blast.pl $results_folder/${filename}/${filename}_precursors.sam $results_folder/${filename}/${filename}_filtered.fa $results_folder/${filename}/${filename}_precursors.fa $var > $results_folder/${filename}/${filename}_precursors.bst 2>>$results_folder/${filename}/${filename}_err
## fi
if [ $bt2tag == "false" ]; then
    echo "Build precursors bowtie index:" >> $results_folder/${filename}/script_err
    echo "bowtie-build -f $results_folder/${filename}/${filename}_precursors.fa $results_folder/${filename}/index/${filename}_precursors" >> $results_folder/${filename}/script_err
	bowtie-build -f $results_folder/${filename}/${filename}_precursors.fa \
                 $results_folder/${filename}/index/${filename}_precursors \
                 >> $results_folder/${filename}/script_log \
                 2>> $results_folder/${filename}/script_err
    echo "Mapping filtered signatures (multiple mapping <= $len) to the precursors:" >> $results_folder/${filename}/script_err
    echo "bowtie -a -v $var -x $results_folder/${filename}/index/${filename}_precursors -f $results_folder/${filename}/${filename}_filtered.fa -p $thread" >> $results_folder/${filename}/script_err
	bowtie -a -v $var \
           -x $results_folder/${filename}/index/${filename}_precursors \
           -f $results_folder/${filename}/${filename}_filtered.fa \
           -p $thread \
           > $results_folder/${filename}/${filename}_precursors.aln \
           2>> $results_folder/${filename}/script_err
	perl $mirdp/convert_bowtie_to_blast.pl $results_folder/${filename}/${filename}_precursors.aln \
         $results_folder/${filename}/${filename}_filtered.fa $results_folder/${filename}/${filename}_precursors.fa > $results_folder/${filename}/${filename}_precursors.bst 2>>$results_folder/${filename}/${filename}_err
fi

sort +3 -25 $results_folder/${filename}/${filename}_precursors.bst > $results_folder/${filename}/${filename}_signatures 2>>$results_folder/${filename}/${filename}_err
echo "finish assigning candidate signatures" $(date +[\ %F\ %X\ ]) >> $results_folder/${filename}/progress_log


#miRDP core algorithm
perl $mirdp/mod-miRDP.pl $results_folder/${filename}/${filename}_signatures $results_folder/${filename}/${filename}_structures > $results_folder/${filename}/${filename}_predictions 2>>$results_folder/${filename}/${filename}_err
echo "finish core algorithm" $(date +[\ %F\ %X\ ]) >> $results_folder/${filename}/progress_log


#calculate chr_length
perl -e '%hash=();@out=();   $name=<>; chomp($name);s/\s+$//; $seq=();  while(<>)  {chomp;  s/\s+$//;  if(m/^[A-Za-z\*]+$/){$seq.=$_;}   elsif(m/^>/){$len=length($seq);$name=~/^>(\S+)\s*/;print "$1\t$len\n"; $name=$_; $seq="";} elsif(m/^$/){next;}  else{print "err! $_\n";} } $len=length($seq);$name=~/^>(\S+)\s*/;print "$1\t$len\n";' $genome > $results_folder/${filename}/chr_length

#remove redundant and filter by plant-specific criteria
perl $mirdp/mod-rm_redundant_meet_plant.pl $results_folder/${filename}/chr_length $results_folder/${filename}/${filename}_precursors.fa $results_folder/${filename}/${filename}_predictions $results_folder/${filename}/${filename}.total_reads $results_folder/${filename}/${filename}_nr_prediction $results_folder/${filename}/${filename}_filter_P_prediction 2>>$results_folder/${filename}/${filename}_err
echo "finish filtering plant criteria" $(date +[\ %F\ %X\ ]) >> $results_folder/${filename}/progress_log


#mkdir $results_folder/${filename}/tmp_files
#mv $results_folder/${filename}/* $results_folder/${filename}/tmp_files 2> /dev/null
#cp $results_folder/${filename}/tmp_files/${filename}_filter_P_prediction $results_folder/${filename}

perl -e ' while(<>){chomp;@tmp=split "\t",$_; $tmp[4]=~m/^(\d+)\.\.(\d+)$/;$mature_beg=$1; $tmp[5]=~m/^(\d+)\.\.(\d+)$/;$pre_beg=$1; if($tmp[1] eq "+"){$pre_beg-=2;$mature_beg-=2;} elsif($tmp[1] eq "-"){$pre_beg+=1;$mature_beg+=1;} $mature_len=length($tmp[6]);$pre_len=length($tmp[7]); $mature_end=$mature_bed+$mature_len;$pre_end=$pre_beg+$pre_len; $sign=""; if(($mature_beg==$pre_beg && $tmp[1] eq "+")||($mature_end==$pre_end && $tmp[1] eq "-") ){$sign="5";}elsif(($mature_beg==$pre_beg && $tmp[1] eq "-")||($mature_end==$pre_end && $tmp[1] eq "+")){$sign="3";}else{$sign="A";}  print "$tmp[0]\t$pre_beg\t$pre_end\t$tmp[3]\t$sign\t$tmp[1]\n";} ' $results_folder/${filename}/${filename}_filter_P_prediction > $results_folder/${filename}/${filename}_filter_P_prediction.bed 2>>$results_folder/${filename}/${filename}_err

rm $formattedInput_filtered $formattedInput

## Reformat the output to include 5p and 3p information
## column names of output are:
## LociName, GenomicCoordinates, 5p mature seq, 5p readsID, 3p mature seq, 3p readsID, precursor seq
awk 'NR==FNR{a[$4]=$1":"$6":"$2}NR>FNR{for(i=1;i<=NF;i++){split($i, tmp, "  \t"); if(length(tmp)==2){k[tmp[1]]=tmp[2]}}; if(k["pri_id"] in a){print k["pri_id"]"\t"a[k["pri_id"]]"\t"k["mature_seq"]"\t"k["star_seq"]"\t"k["pre_seq"]}}' FS="\t" RS="\n" $results_folder/${filename}/${filename}_filter_P_prediction RS="\n\n" FS="\n" $results_folder/${filename}/${filename}_predictions | # Extract MIR locus information from miRDeep result
    awk '{k=index($5,$3); if(k<length($5)/2){print $1"\t"$2"\t"$3"\t"$4"\t"$5}else{print $1"\t"$2"\t"$4"\t"$3"\t"$5}}' | # identify 5p and 3p
    awk 'NR==FNR{a[$2]=$1}NR>FNR{if($3 in a){p5=a[$3]}else{p5="NA"}; if($4 in a){p3=a[$4]}else{p3="NA"}; print $1"\t"$2"\t"$3"\t"p5"\t"$4"\t"p3"\t"$5}' <(fasta_formatter -t -i $results_folder/${filename}/${filename}.fa) - > $results_folder/${filename}/${filename}_filter_P_prediction_reformatted

## If refFile provided, miRNA will be compared with reference set, and be named accordingly
## The refFile should contain 6 columns: 5p_seq, 3p_seq, 5p_name, 3p_name, precursor_name, precursor_seq
## all known miRNA should be sorted and uniqufied by 5p_seq, 3p_seq (bedtools groupby)

## example for prepare refFile from PmiREN2.0
## soybean known miRs
## cut -f 1,13,15,18,19,22 Glycine_max_basicInfo.txt |
##     awk '
##     BEGIN{OFS="\t"}
##     NR==1{print}
##     NR>1{k=length($2); sub("*","", $5); if(index($2,$4)<=k/2){$3=$3"-5p";$5=$5"-3p"}else{$3=$3"-3p";$5=$5"-5p"}; print}
##     ' |
##     awk 'NR>1{if($3~/5p$/){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}else{print $1"\t"$2"\t"$5"\t"$6"\t"$3"\t"$4}}' |
##     awk '{print $3"\t"$4"\t"$5"\t"$6"\t"$1"\t"$2}' |
##     sort -k 2,2 -k 4,4 -k 1,1 -k 3,3 |
##     bedtools groupby -g 2,4 -c 1,3,5,6 -o collapse > mature_precursor_seq.tsv

awk '
BEGIN{FS="\t"; OFS="\t"}
NR==FNR{
    a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6
}NR>FNR{
    fiveMatureName="";
    threeMatureName="";
    preName="";
    if($3"\t"$5 in a){
        split(a[$3"\t"$5], k, "\t");
        fiveMatureName=k[1];
        threeMatureName=k[2];
        preName=k[3];
        print "Identified\t"$1"\t"$2"\t"fiveMatureName"\t"$3"\t"$4"\t"threeMatureName"\t"$5"\t"$6"\t"preName"\t"$7
    }else{
        for(i in a){
            split(i, p, "\t");
            split(a[i], q, "\t");
            fiveMatureRef=p[1];threeMatureRef=p[2];
            split(q[4], preSeqRef, ",")
            split(q[3], preNameRef, ",")
            if(fiveMatureRef==$3){
                fiveMatureName=q[1];
                for(m in preSeqRef){
                    if(index(preSeqRef[m], $5)){
                        threeMatureName="["q[2]"]-iso"
                        preName=preName","preNameRef[m]
                    }
                }
                break
            }else if(threeMaturename==$5){
                threeMatureName=q[2];
                for(m in preSeqRef){
                    if(index(preSeqRef[m], $3)){
                        fiveMatureName="["q[1]"]-iso"
                        preName=preName","preNameRef[m]
                    }
                }
                break
            }else{
                for(m in preSeqRef){
                    if(index(preSeqRef[m], $3) && index(preSeqRef[m], $5)){
                        fiveMatureName="["q[1]"]-iso"
                        threeMatureName="["q[2]"]-iso"
                        preName=preName","preNameRef[m]
                    }
                }
            }
        }
        if(fiveMatureName==""){fiveMatureName="Novel"}
        if(threeMatureName==""){threeMatureName="Novel"}
        if(preName==""){preName="Novel"}else{preName=substr(preName,2)}
        print "Iso/Novel\t"$1"\t"$2"\t"fiveMatureName"\t"$3"\t"$4"\t"threeMatureName"\t"$5"\t"$6"\t"preName"\t"$7
    }
}
' $refFile $results_folder/${filename}/${filename}_filter_P_prediction_reformatted > $results_folder/${filename}/${filename}_filter_P_prediction_reformatted_withFamilyInfo.tsv
echo "finish reformat and family annotation" $(date +[\ %F\ %X\ ]) >> $results_folder/${filename}/progress_log

expr_c1=$(mktemp -p $results_folder)
awk '{if($6!="NA"){print $2"-5p\t"$4"\t"$5}; if($9!="NA"){print $2"-3p\t"$7"\t"$8}}' $results_folder/${filename}/${filename}_filter_P_prediction_reformatted_withFamilyInfo.tsv > $expr_c1

function process_each_sample {

    local inputFile=$1
    local inputFormat=$2
    local inputProcessed=$(mktemp -p $results_folder)
    local ncRNA_out=$(mktemp -p $results_folder)
    local bowtieFormatTag=""
    local inputFileName=$(basename $inputFile)
    local sampleName=${inputFileName%%.[fq|fa]*}

    if [[ $inputFile =~ .(fa|fasta).gz$ ]] || [[ $inputFile =~ .(fa|fasta)$ ]]; then
        bowtieFormatTag="-f"
    elif [[ $inputFile =~ .(fq|fastq).gz$ ]] || [[ $inputFile =~ .(fq|fastq)$ ]]; then
        bowtieFormatTag="-q"
    elif [[ ! -z $inputFormat ]]; then
        bowtieFormatTag=$inputFormat
    else
        echo "The input format should be .fa, .fasta, .fq or .fastq, or formatTag should be provided."
        exit 1
    fi
    ## map to rfam_ncRNA, remove reads with hits
    bowtie -v 0 \
           -x $mirdp/index/rfam_index \
           -p $thread \
           $bowtieFormatTag $inputFile |
            awk '{print $1}' |
            sort -u > $ncRNA_out
    bowtie -v $var $large \
           -m $len \
           -x $bowtie_index \
           -p $thread \
           -S \
           $bowtieFormatTag $inputFile |
        awk '
            NR==FNR{a[$1]}
            NR>FNR{if($1~/^@/){print}else{if(!($1 in a)){print}}}
            ' $ncRNA_out - |
        samtools view -@ $thread -b > $results_folder/${filename}/${sampleName}.bam

    if [[ $inputFile =~ .gz$ ]]; then
        zcat $inputFile | fastx_collapser -Q33  > $inputProcessed
    else
        fastx_collapser -Q33 -i $inputFile > $inputProcessed
    fi
    local lib_size=$(samtools view $results_folder/${filename}/${sampleName}.bam | awk '$3!="*"{print $1}' | sort -u |wc -l)
    awk 'NR==FNR{split($1, tmp, "-"); a[$2]=tmp[2]}NR>FNR{if($3 in a){print $0"\t"a[$3]/'$lib_size'*1000000}else{print $0"\t0"}}' <(fasta_formatter -t -i $inputProcessed) $expr_c1 > $results_folder/${filename}/${sampleName}.expr.tsv
    rm $ncRNA_out $inputProcessed
}

if [[ ! -z $input ]] || [[ -f $batch ]]; then
    if [[ $input =~ .gz$ && tag == "f" ]]; then
        ## assume formatted input provided, do not estimate expression
        true # PASS
    elif [[ $input =~ .(fa|fasta)$ && $tag == "f" ]]; then
        ## assume formatted input provided, do not estimate expression
        true # PASS
    elif [[ $input =~ .(fq|fastq).gz$ && $tag == "q" ]]; then
        ## Assume single library fastq
        if [[ $trim == "true" ]]; then
            process_each_sample $inputTrimmed "-f"
        else
            process_each_sample $input
        fi
    elif [[ $input =~ .(fq|fastq)$ && $tag == "q" ]]; then
        ## Assume single library fastq
        if [[ $trim == "true" ]]; then
            process_each_sample $inputTrimmed "-f"
        else
            process_each_sample $input
        fi
    elif [[ -f $batch ]]; then
        ## Multiple libraries fasta/fastq(.gz)
        if [[ $trim == "true" ]]; then
            for i in $(cat $batchTrimmed); do process_each_sample $i "-f"; done
        else
            for i in $(cat $batch); do process_each_sample $i; done
        fi
    else
        echo "Either '--input' or '--batch' option is not valid."
        exit 1
    fi
fi

## Remove trimmed files
if [[ ! -z $batchTrimmed ]]; then
    awk 'BEGIN{a=""}{a=a" "$1}END{print "rm "a}' $batchTrimmed |
        bash
    rm $batchTrimmed
else
    rm $inputTrimmed
fi

## Merge Expr table
tmpExpr1=$(mktemp -p $results_folder)
tmpExpr2=$(mktemp -p $results_folder)
awk 'BEGIN{print "id\tfamily\tmature_seq"}{print}' $expr_c1 > $tmpExpr1
for i in $results_folder/${filename}/*.expr.tsv;
do
    sampleID=$(basename "$i" .expr.tsv)
    awk '
    NR==FNR{a[$1"\t"$2"\t"$3]=$4}
    NR>FNR&&FNR==1{print $0"\t'$sampleID'"}
    NR>FNR&&FNR>1{print $0"\t"a[$1"\t"$2"\t"$3]}
    ' $i $tmpExpr1 > $tmpExpr2;
    mv $tmpExpr2 $tmpExpr1;
done
mv $tmpExpr1 $results_folder/${filename}/combined.expr.tsv
rm $expr_c1

echo "pipeline finished" $(date +[\ %F\ %X\ ]) >> $results_folder/${filename}/progress_log


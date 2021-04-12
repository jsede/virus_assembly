#!/bin/bash

Version=virus_assembly_V1.0

USAGE="
Pipeline: NGS pipeline for viral assembly.
usage: $(basename "$0" .sh) [-h -v -p -q] (-i dir -m value -t value )
(-s string) 
with:
    -h  Show help text
    -v  Version of the pipeline
    -n  Name of RUN.
    -i  Input directory
    -s  Viral species [HIV, RRV, HMPV]
    -c  Perform clipping of primers
    -q  Perform quality check using fastQC
    -m  Memory
    -t  Number of threads
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hvi:o:s:pqm:t:n:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
	v) echo "$Version"; exit 1;;
    i) IN_DIR=$OPTARG;;
    s) VIRUS=$OPTARG;;
    c) CLIP_FLAG=YES;;
    q) QC_FLAG=YES;;
    m) MEMORY=$OPTARG;;
    n) RUN=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done




# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${IN_DIR+x} ]; then echo "-i $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${VIRUS+x} ]; then echo "-s $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RUN+x} ]; then echo "-n is missing. Defaulting to OUTPUT"; echo ""; RUN="OUTPUT"; fi; 
if [ -z ${MEMORY+x} ]; then echo "-m is missing. Defaulting to ?"; echo ""; MEMORY=2; fi;
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 2 threads."; THREADS=2; fi;


### Define values ---------------------------------------------------------

#define project folder

#set database folder
db=$CONDA_PREFIX/virus_assembly/db

#set virus for filtering
blast="$VIRUS"_refs

#set compute resources where m=GB t=CPU
m=$MEMORY
t=$THREADS

#change into project folder
cd $IN_DIR

#quality trim reads with bbduk/37.99
while read i; do
bbduk.sh -Xmx"$MEMORY"g threads="$THREADS" in=$IN_DIR/1_reads/"$i"_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_2.fastq.gz out=stdout.fq ref=phix k=31 hdist=1 | bbduk.sh -Xmx"$MEMORY"g threads="$THREADS" interleaved=true in=stdin.fq out=stdout.fq ref=adapters ktrim=r k=17 mink=3 | bbduk.sh -Xmx"$MEMORY"g threads="$THREADS" interleaved=true in=stdin.fq out=$IN_DIR/1_reads/"$i"_trim_1.fastq.gz out2=$IN_DIR/1_reads/"$i"_trim_2.fastq.gz minlen=50 qtrim=rl trimq=20 entropy=0.7
done < $IN_DIR/IDs.list

#normalise read coverage for denovo assembly
while read i; do
bbnorm.sh -Xmx"$MEMORY"g threads="$THREADS" in=$IN_DIR/1_reads/"$i"_trim_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_2.fastq.gz out=$IN_DIR/1_reads/"$i"_trim_norm_1.fastq.gz out2=$IN_DIR/1_reads/"$i"_trim_norm_2.fastq.gz target=100;
done < $IN_DIR/IDs.list

#QC check
if [ "${QC_FLAG}" == "YES" ]; then 
while read i; do
fastqc -t "$THREADS" $IN_DIR/1_reads/"$i"_*.fastq.gz
done < $IN_DIR/IDs.list
multiqc $IN_DIR/1_reads/*_fastqc.zip --ignore *_norm* -o $IN_DIR/1_reads/multiqc_fastQC_report -n multiqc_fastQC_report
rm $IN_DIR/1_reads/*_fastqc.zip $IN_DIR/1_reads/*_fastqc.html; fi;

#CLIPPING
if [ "${CLIP_FLAG}" == "YES" ]; then 
while read i; do
cutadapt -j "$THREADS" -a file:$db/primers/"$VIRUS".primers.clips.A.fa -A file:$db/primers/"$VIRUS".primers.clips.A.fa -o $IN_DIR/1_reads/"$i"_trim_temp1_1.fastq.gz -p $IN_DIR/1_reads/"$i"_trim_temp1_2.fastq.gz $IN_DIR/1_reads/"$i"_trim_1.fastq.gz $IN_DIR/1_reads/"$i"_trim_2.fastq.gz
cutadapt -j "$THREADS" -g file:$db/primers/"$VIRUS".primers.clips.B.fa -G file:$db/primers/"$VIRUS".primers.clips.B.fa -o $IN_DIR/1_reads/"$i"_trim_temp2_1.fastq.gz -p $IN_DIR/1_reads/"$i"_trim_temp2_2.fastq.gz $IN_DIR/1_reads/"$i"_trim_temp1_1.fastq.gz $IN_DIR/1_reads/"$i"_trim_temp1_2.fastq.gz
bbduk.sh -Xmx"$MEMORY"g threads="$THREADS" in=$IN_DIR/1_reads/"$i"_trim_temp2_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_temp2_2.fastq.gz out=$IN_DIR/1_reads/"$i"_trim_clip_1.fastq.gz out2=$IN_DIR/1_reads/"$i"_trim_clip_2.fastq.gz minlen=50
rm $IN_DIR/1_reads/"$i"_trim_temp*
done < $IN_DIR/IDs.list; fi;

#create folder for reference mapping
mkdir $IN_DIR/2_ref_map

#map trimmed read to reference genomes
while read i; do
    if [ "${CLIP_FLAG}" == "YES" ]; 
        then bbmap.sh -Xmx"$MEMORY"g threads="$THREADS" maxindel=500 in=$IN_DIR/1_reads/"$i"_trim_clip_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_clip_2.fastq.gz outm=$IN_DIR/2_ref_map/"$i".ref_mapped.bam ref=$db/fasta/"$VIRUS"_refs.fa;
else bbmap.sh -Xmx"$MEMORY"g threads="$THREADS" maxindel=500 in=$IN_DIR/1_reads/"$i"_trim_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_2.fastq.gz outm=$IN_DIR/2_ref_map/"$i".ref_mapped.bam ref=$db/fasta/"$VIRUS"_refs.fa;
fi;
samtools sort -@ "$THREADS" -o $IN_DIR/2_ref_map/"$i".ref_mapped.sorted.bam $IN_DIR/2_ref_map/"$i".ref_mapped.bam
rm $IN_DIR/2_ref_map/"$i".ref_mapped.bam
done < $IN_DIR/IDs.list

#generate coverage maps
cd $IN_DIR/2_ref_map
while read i; do
qualimap bamqc -bam "$i".ref_mapped.sorted.bam  -outformat PDF -outfile "$i".ref_mapped.coverage.pdf --java-mem-size="$MEMORY"G
done < $IN_DIR/IDs.list
cd $IN_DIR

#cleanup files
while read i; do
mv $IN_DIR/2_ref_map/"$i".ref_mapped.sorted_stats/"$i".ref_mapped.coverage.pdf $IN_DIR/2_ref_map/
rm -r $IN_DIR/2_ref_map/"$i".ref_mapped.sorted_stats
done < $IN_DIR/IDs.list

#create folder for de novo assembly
mkdir $IN_DIR/3_contigs

#denovo assemble trimmed reads with megahit and rename output contigs with library name
while read i; do
megahit -t "$THREADS" -1 $IN_DIR/1_reads/"$i"_trim_norm_1.fastq.gz -2 $IN_DIR/1_reads/"$i"_trim_norm_2.fastq.gz -o $IN_DIR/3_contigs/megahit_"$i" --out-prefix "$i" --min-contig-len 500;
sed "s/>k/>Draft_"$i"_k/g" $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.fa > $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.fa;
done < $IN_DIR/IDs.list

#export denovo assembly contig stats
while read i; do
grep ">" $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.fa | cut -d ">" -f2 | tr ' ' '\t' | sed 's/flag\=//g' | sed 's/multi\=//g' | sed 's/len\=//g' > $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.stats
done < $IN_DIR/IDs.list

#extract high coverage denovo assembly contigs
while read i; do
awk -F$'\t' '{OFS=FS}{if ($3>50) print $1}' $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.stats > $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.hicov.list
seqtk subseq $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.fa $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.hicov.list > $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.hicov.fa
done < $IN_DIR/IDs.list

#create folder for blasting
mkdir $IN_DIR/4_filter

#combine denovo contigs into single file
cat $IN_DIR/3_contigs/megahit_*/*.contigs.megahit.hicov.fa > $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.fasta

#blast to local database
export BLASTDB=$db/blast
blastn -query $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.fasta -db "$db"/blast/"$VIRUS"_refs -evalue 1E-10 -num_threads 1 -out $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$VIRUS".txt -outfmt "6 qseqid qlen stitle sstart send pident length evalue sstrand"

#get top virus blast results for each contig
awk -F$'\t' '!seen[$1]++' $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$VIRUS".txt > $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$VIRUS".top.txt

#take column containing contig names
cut -f1 $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$VIRUS".top.txt > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".list 

#retrieve sequences
seqtk subseq $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.fasta $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".list > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".fa

#align to reference
mafft --thread "$THREADS" --reorder --adjustdirection --maxiterate 10 --add $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".fa $db/fasta/"$VIRUS"_refs.fa > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_aligned.fa

#unalign draft genome alignment
sed '/^>/! s/\-//g' $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_aligned.fa | sed 's/_R_//g' | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_unaligned.fa

#retrieve draft denovo genomes
seqtk subseq $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_unaligned.fa $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".list > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_unaligned.reoriented.fa

#create folder for blasting
mkdir $IN_DIR/5_remap

#retrieve individual draft sequences for mapping
while read i; do
grep -A1 ">Draft_${i}_" $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_unaligned.reoriented.fa > $IN_DIR/5_remap/"$i".draft.fa
done < $IN_DIR/IDs.list

#map trimmed (clipped) reads to draft genomes
while read i; do
    if [ "${CLIP_FLAG}" == "YES" ]; then 
        bbmap.sh -Xmx"$MEMORY"g threads="$THREADS" maxindel=200 minid=0.98 in=$IN_DIR/1_reads/"$i"_trim_clip_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_clip_2.fastq.gz outm=$IN_DIR/5_remap/"$i".remapped.bam ref=$IN_DIR/5_remap/"$i".draft.fa;
else bbmap.sh -Xmx"$MEMORY"g threads="$THREADS" maxindel=200 minid=0.98 in=$IN_DIR/1_reads/"$i"_trim_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_2.fastq.gz outm=$IN_DIR/5_remap/"$i".remapped.bam ref=$IN_DIR/5_remap/"$i".draft.fa
fi;
samtools sort -@ "$THREADS" -o $IN_DIR/5_remap/"$i".remapped.sorted.bam $IN_DIR/5_remap/"$i".remapped.bam
rm $IN_DIR/5_remap/"$i".remapped.bam
done < $IN_DIR/IDs.list

# Used scripts


```bash
#1.remove_small.pl
#2.header.fasta.numbers.pl
#3.hitter.py
#4.hitter_table.py
#5.m5nr_annotation.py
#6.hitter_na.py 
#7.promer_deid_v9gmv.py
```

# Host filtering sequences 


```bash
#Create host genome index database with Bowtie2

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
bowtie2-build tomato_genome.fna tomato_genome
```


```bash
#Map tomato rhizosphere metagenomes against host (tomato) genome with Bowtie2

#Align against reference genome
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

SEQS=host_filtering


for FAA in `ls *.fastq.gz | perl -pe 's/\_.*//g' | sort | uniq`
do

echo "#!/bin/bash" >hostmap_$FAA.sh
echo "#$ -cwd" >>hostmap_$FAA.sh
echo "#$ -j y" >>hostmap_$FAA.sh
echo "#$ -S /bin/bash" >>hostmap_$FAA.sh
echo "bowtie2 -x tomato_genome -1 $SEQS/$FAA"_R1.fastq.gz" -2 $SEQS/$FAA"_R2.fastq.gz" -S $FAA.sam"  -p 20 >>hostmap_$FAA.sh
chmod +x *.sh; done
```


```bash
#Convert sam to bam

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

FILES=host_filtering

for FAA in `ls *.sam | sed -e 's/.sam//g'`
do
echo "#!/bin/bash" >samtobam__$FAA.sh
echo "#$ -cwd" >>samtobam__$FAA.sh
echo "#$ -j y" >>samtobam__$FAA.sh
echo "#$ -S /bin/bash" >>samtobam__$FAA.sh
echo "#!/bin/bash" >>samtobam__$FAA.sh
echo "samtools view -bS $FILES/$FAA".sam "> $FAA.bam" >>samtobam__$FAA.sh
chmod +x *.sh; done
```


```bash
#Get unmaped paired sequences against tomato genome

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

FILES=host_filtering

for FAA in `ls *.bam | sed -e 's/.bam//g'`
do

echo "#!/bin/bash" >getunmaped_$FAA.sh
echo "#$ -cwd" >>getunmaped_$FAA.sh
echo "#$ -j y" >>getunmaped_$FAA.sh
echo "#$ -S /bin/bash" >>getunmaped_$FAA.sh
echo "#!/bin/bash" >>getunmaped_$FAA.sh

echo "samtools view -b -f 12 -F 256  $FILES/$FAA".bam  "> $FAA"_umh_bothends.bam  >>getunmaped_$FAA.sh
chmod +x *.sh; done
```


```bash
#Obtain R1 and R2 paired sequences for host filtered samples
#sort bam file by read name (-n) to have paired reads next to each other as required by bedtools


#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#!/bin/bash

FILES=host_filtering

for FAA in `ls *_umh_bothends.bam | sed -e 's/_umh_bothends.bam//g'`
do

echo "#!/bin/bash" >sortedbam_$FAA.sh
echo "#$ -cwd" >>sortedbam_$FAA.sh
echo "#$ -j y" >>sortedbam_$FAA.sh
echo "#$ -S /bin/bash" >>sortedbam_$FAA.sh
echo "#!/bin/bash" >>sortedbam_$FAA.sh

echo "samtools sort -n $FILES/$FAA"_umh_bothends.bam -o "$FAA"_sorted.bam  >>sortedbam_$FAA.sh
chmod +x *.sh; done
```


```bash
#Get the filtered fastq files against host genome using bedtools

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#!/bin/bash

FILES=host_filtering
BIN=bin


for FAA in `ls *_sorted.bam  | sed -e 's/_sorted.bam//g'`
do
echo "#!/bin/bash" >filth_$FAA.sh
echo "#$ -cwd" >>filth_$FAA.sh
echo "#$ -j y" >>filth_$FAA.sh
echo "#$ -S /bin/bash" >>filth_$FAA.sh

echo "$BIN"/bamToFastq  -i  "$FILES/$FAA"_sorted.bam "-fq $FILES/$FAA"_filtered_R1.fastq "-fq2 $FILES/$FAA"_filtered_R2.fastq >> filth_$FAA.sh
chmod +x *.sh; done
```

# Quality filtering


```bash
#Using Trimmomatic for Quality Filtering and Processing of Paired and Unpaired Sequences

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

SEQS=trimming
BIN=bin/Trimmomatic-0.36

COUNT=0
for FAA in `ls *.fastq| perl -pe 's/\_.*//g' | sort | uniq`
do
echo "#!/bin/bash" >trim_$FAA.sh
echo "#$ -cwd" >>trim_$FAA.sh
echo "#$ -j y" >>trim_$FAA.sh
echo "#$ -S /bin/bash" >>trim_$FAA.sh
echo "java -jar $BIN/trimmomatic-0.36.jar PE -threads 2 -phred33 -trimlog $SEQS/trim.log $SEQS/$FAA"_filtered_R1.fastq" $SEQS/$FAA"_filtered_R2.fastq" $SEQS/$FAA"_paired_R1.fastq.gz" $SEQS/$FAA"_unpaired_R1.fastq.gz "$SEQS/$FAA"_paired_R2.fastq.gz" $SEQS/$FAA"_unpaired_R2.fastq.gz" ILLUMINACLIP:/bin/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" >>trim_$FAA.sh
chmod +x *.sh; done
```

# Hybrid Assembly 


```bash
#Assembly Spades

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

SEQS=hybrid_ass
BIN=bin/SPAdes-3.12.0-Linux/bin

COUNT=0
for FAA in `ls *.fastq.gz | perl -pe 's/\_.*//g' | sort | uniq`
do
echo "#!/bin/bash" >$FAA.scr
echo "#$ -cwd" >>$FAA.scr
echo "#$ -j y" >>$FAA.scr
echo "#$ -S /bin/bash" >>$FAA.scr
echo date >> $FAA.scr
echo  "$BIN"/spades.py --meta -1 "$SEQS/$FAA"_paired_R1.fastq.gz -2 "$SEQS/$FAA"_paired_R2.fastq.gz -o "$SEQS/contigs_spades/contig_$FAA"/ >>$FAA.scr
echo date >>$FAA.scr
chmod +x *.scr; done
```


```bash
#Mapping unassembled sequences from SPAdes to create a second assembly using velvet

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

CONT=hybrid_ass/map_spades
COUNT=0

for FAA in `ls *paired_R1* | sed -e 's/_paired_R1.fastq.gz//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >spadesmap_$FAA.scr
echo "#$ -cwd" >>spadesmap_$FAA.scr
echo "#$ -j y" >>spadesmap_$FAA.scr
echo "#$ -S /bin/bash" >>spadesmap_$FAA.scr

echo bbwrap.sh ref="$CONT"/contigs_spades_"$FAA".fasta in="$CONT/$FAA"_paired_R1.fastq.gz in2="$CONT/$FAA"_paired_R2.fastq.gz build="$COUNT" out="$FAA".sam kfilter=22 subfilter=15 maxindel=80 >>spadesmap_$FAA.scr
chmod +x *.scr; done
```


```bash
#Get Spades unmapped reads

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

SEQS=hybrid_ass/map_spades

for FAA in `ls *.sam | sed -e 's/.sam//g'`
do
echo "#!/bin/bash" >get_un_spa_$FAA.scr
echo "#$ -cwd" >>get_un_spa_$FAA.scr
echo "#$ -j y" >>get_un_spa_$FAA.scr
echo "#$ -S /bin/bash" >>get_un_spa_$FAA.scr

echo samtools view -u -f 12 -F 256 "$SEQS/$FAA".sam "|" samtools bam2fq - ">" "$FAA"_ump.fastq  >> get_un_spa_$FAA.scr
chmod +x *.scr; done
```


```bash
#Separation of unmapped reads into forward and reverse

for FAA in `ls *_ump.fastq | sed -e 's/_ump.fastq//g'`
do
grep '^@.*/1$' -A 3 --no-group-separator "$FAA"_ump.fastq > "$FAA"_um_r1.fastq
grep '^@.*/2$' -A 3 --no-group-separator "$FAA"_ump.fastq > "$FAA"_um_r2.fastq; done

```


```bash
#Assembly Velvet
#Assembling reads that did not map to Spades contigs


#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

SEQS=hybrid_ass/velvet
BIN=bin/velvet

COUNT=0
for FAA in `ls *_um_r1.fastq | perl -pe 's/\_um_r1.fastq//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >velvet_$FAA.scr
echo "#$ -cwd" >>velvet_$FAA.scr
echo "#$ -j y" >>velvet_$FAA.scr
echo "#$ -S /bin/bash" >>velvet_$FAA.scr


echo "$BIN"/velveth "$SEQS/$FAA"_assemblyVelvet 31 -shortPaired -fastq "$FAA"_um_r1.fastq "$FAA"_um_r2.fastq >>velvet_$FAA.scr
echo "$BIN"/velvetg "$SEQS/$FAA"_assemblyVelvet -exp_cov 2 -ins_length 350 >> velvet_$FAA.scr
chmod +x *.scr; done
```


```bash
# Merge of Spades and Velvet contigs 
for FAA in `ls contigs_spades* | sed -e 's/contigs_spades_//g' |  sed -e 's/.fasta//g'`
do
cat /hybrid_ass/contigs_spades/all_contigs_spades/contigs_spades_"$FAA".fasta /hybrid_ass/velvet/all_contigs_velvet/contigs_velvet_"$FAA".fasta> assh_merge/"$FAA"_hb_contigs.fasta;done

```


```bash
#Remove reads shorter than 100pb
for FAA in `ls *_hb_contigs.fasta| sed -e 's/_hb_contigs.fasta//g'`; do perl remove_small.pl 100 "$FAA"_hb_contigs.fasta> "$FAA"_hb100_contigs.fasta; done

```


```bash
# Index hydrid contigs according to sample name
for i in `ls *_hb100_contigs.fasta| sed -e 's/_hb100_contigs.fasta//g'`; do perl header.fasta.numbers.pl "$i" "$i"_hb100_contigs.fasta; done

```

# ORFs prediction


```bash
#Predicting ORFs from contigs usingprodigal

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

SEQS=orf_predict
COUNT=0
for FAA in `ls *.numbered.fas | sed -e 's/\_hb100_contigs.fasta.numbered.fas//g' | sort | uniq`
do
echo "#!/bin/bash" >prodigal_"$FAA"_.scr
echo "#$ -cwd" >>prodigal_"$FAA"_.scr
echo "#$ -j y" >>prodigal_"$FAA"_.scr
echo "#$ -S /bin/bash" >>prodigal_"$FAA"_.scr
echo /bin/Prodigal/prodigal -i "$SEQS/$FAA"_hb100_contigs.fasta.numbered.fas -o  "$SEQS/$FAA"_prodigal.out -d "$SEQS/$FAA".fna -a  "$SEQS/$FAA".faa -p meta >>prodigal_"$FAA"_.scr
chmod +x *.scr; done
```

# Predicted proteins annotation against M5nr database


```bash
#Split predicted proteins file to parallelize the annotation 
partefasta 72500 hydro2.faa 
partefasta 17250 hydro3.faa 
partefasta 40600 hydro4.faa 
partefasta 30000 soil1.faa
partefasta 32000 soil2.faa 
partefasta 24500 soil3.faa 
partefasta 31500 soil4.faa 
partefasta 26500 soil5.faa 
```


```bash
#Aligning predicted proteins against the M5nr database using diamond 

#!/bin/bash
#bash nombre_shipt.sh <nombre-del-trabajo>
SEQS=anotation_m5nr
BIN=bin
DB=databases/m5NRv2020
COUNT=0

for FAA in `ls *.fas`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo  "$BIN/diamond blastp -d $DB/m5nr_LDAdic20.dmnd -q $SEQS/$FAA" -f 6 -e 1e-10 -k 10 -p 1 -o $SEQS/"$FAA".bout  >>$*.$COUNT.scr
chmod +x *.scr; done

```


```bash
#Concatenate all outputs belonged to the same sample
for N in `ls *.faa | sed -e 's/\.faa//g' | sort | uniq`;do cat "$N".faa*.bout >"$N"_m5nr.bout;done
```


```bash
#Get the best hits from diamond output tables based on bitscore values
for FAA in `ls *m5nr.bout | sed -e 's/_m5nr.bout//g'`; do cat "$FAA"_m5nr.bout |  perl -pe '$name_col=0; $score_col=11; while(<>) { s/\r?\n//; @F=split /\t/, $_; ($n, $s) = @F[$name_col, $score_col]; if (! exists($max{$n})) { push @names, $n }; if (! exists($max{$n}) || $s > $max{$n}) { $max{$n} = $s; $best{$n} = () }; if ($s == $max{$n}) { $best{$n} .= "$_\n" }; } for $n (@names) { print $best{$n} } ' > best;  perl -e ' $column=0; $unique=0; while(<>) { s/\r?\n//; @F=split /\t/, $_; if (! ($save{$F[$column]}++)) { print "$_\n"; $unique++ } } ' best > "$FAA"_best_uniq; rm best; done

#Filter in order to get only 60% identity
for N in `ls *_best_uniq| sed -e 's/\_best_uniq//g'`; do awk -F "\t" '{if ($3>60) print $0}' "$N"_best_uniq > "$N"_pident_m5nr.bout;done

#Simplify output of best hits files
for i in `ls *_pident_best_uniq| sed -e 's/\_pident_best_uniq//g'`; do awk '{print $1"\t"$3"\t"$11"\t"$12"\t"$2}' "$i"_pident_best_uniq > "$i"_best.simple.tsv; done

```

# Processing metagenomes subsamples


```bash
#Extracting 1,000,000 random sequences from each metagenome

# Define the sample names and random seeds
samples=("hydro2" "hydro3" "hydro4" "soil1" "soil2" "soil3" "soil4" "soil5")
seeds=(123 354 24 588 69 157 47 963 454 52 498 862 751 269 342 419 958 215 1058 4569 65 549 784 4495)

# Loop through each sample and seed set 
for sample in "${samples[@]}"; do
  for i in {1..3}; do
    seed=${seeds[$((3 * (${samples[@]/$sample//} / (4 * ${#samples[@]})) + $i - 1))]}
    seqtk sample -s $seed ${sample}_paired_R1.fastq 1000000 > ${sample}_ss${i}_R1.fastq
    seqtk sample -s $seed ${sample}_paired_R2.fastq 1000000 > ${sample}_ss${i}_R2.fastq
  done
done
```


```bash
#Covert fastq to fasta 
for N in `ls *_ss*_R1.fastq| sed -e 's/_R1.fastq//g'`; do sed -n '1~4s/^@/>/p;2~4p' "$N"_R1.fastq > "$N"_R1.fasta; done
for N in `ls *_ss*_R2.fastq| sed -e 's/_R2.fastq//g'`; do sed -n '1~4s/^@/>/p;2~4p' "$N"_R2.fastq > "$N"_R2.fasta; done

#cat files belonging to the same sample 
for FAA in `ls *_R1.fasta | sed -e 's/_R1.fasta//g'`; do cat "$FAA"*.fasta > "$FAA"_rn.fasta;done
```

# Mapping reads against predicted proteins to estimate abundance


```bash
#Creating databases for each ORF prediction file

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

SEQS=subsamples

for FAA in `ls *.fna | sed -e 's/.fna//g'`
do
echo "#!/bin/bash" > bbuild_"$FAA".scr
echo "#$ -cwd" >>bbuild_"$FAA".scr
echo "#$ -j y" >>bbuild_"$FAA".scr
echo "#$ -S /bin/bash" >>bbuild_"$FAA".scr
echo bowtie2-build "$SEQS/$FAA".fna "$SEQS/$FAA">>bbuild_"$FAA".scr
chmod +x *.scr; done
```


```bash
#Mapping reads against samples and subsamples ORFs

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

SEQS=subsamples

for FAA in `ls {hydro*,soil*}.fasta* 2>/dev/null | sed -e 's/\.fasta.*//g'`
do
    echo "#!/bin/bash" > bowtie_"$FAA".scr
    echo "#$ -cwd" >> bowtie_"$FAA".scr
    echo "#$ -j y" >> bowtie_"$FAA".scr
    echo "#$ -S /bin/bash" >> bowtie_"$FAA".scr

    echo "bowtie2 -f -x $FAA -U \"$SEQS/$FAA\"ss*.fasta -S \"$SEQS/$FAA".sam\" --quiet -p 20 --very-sensitive" >> bowtie_"$FAA".scr
    chmod +x bowtie_"$FAA".scr
done
```


```bash
#Obtain quality alignments to reference sequences and their frequency
for N in `ls *.sam | sed -e 's/.sam//g'`; do grep -v '^@' "$N".sam | awk '{if($5 == "42") print $3}' | sort | uniq -c > "$N".hits;done
#Make table
for FAA in `ls *.fna | sed -e 's/.fna//g'`; do grep ">" "$FAA".fna | sed 's/#/\t/g;s/>//g' | cut -f1 | sed 's/ /\t/g' | awk 'BEGIN{i=0} /.*/{printf "%d\t% s\n",i,$0; i++}' | cut -f1,2 > "$FAA".otu; done

#list_hits contains names of output files from mapping with bowtie2
#list contains names of samples and subsamples

```


```bash
#Using hitter.py script to obtain abundances of each protein (hits to M5nr) by parsing diamond output tables and counting hits per sample

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

OUT=abundances/hitter
BIN=bin
ANT=anotation_m5nr
ALIG=subsamples

for HIT in `cat list_hits | sed -e 's/.hits//g'`
do
for FAA in `cat lista`
do
echo "#!/bin/bash" >hitter_"$FAA"_"$HIT".scr
echo "#$ -cwd" >>hitter_"$FAA"_"$HIT".scr
echo "#$ -j y" >>hitter_"$FAA"_"$HIT".scr
echo "#$ -S /bin/bash" >>hitter_"$FAA"_"$HIT".scr
echo python3 "$BIN"/hitter.py "$ANT/$FAA"_best.simple.tsv "$ALIG/$HIT".hits "$ALIG/$HIT".otu "$HIT"  >>hitter_"$FAA"_"$HIT".scr
done
chmod +x *.scr; done


#list_hits contains names of output files from mapping with bowtie2
#list contains names of samples and subsamples
```


```bash
#Joining hits tables with the script hitter_table.py

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

python3 hitter_table.py lista_hout hs_subsamples_global

#lista_hout contains names of output files hitter.py
#The output of hitter_table.py contains a table named hitter_table.py with the abundance of each M5nr identifier in each sample
```

# Generate tables of predicted proteins annotation and counts

## Refseq


```bash
#Obtain  M5nr identifiers for annotation against Refseq adatabase
cut -f1 hs_subsamples_global.tsv | sed '1d' >ids_md5_subsamples.txt

#split MD5 identifiers
split -l 198755 ids_md5_subsamples.txt ids_md5_subsamples_split_
```


```bash
#Obtain unique M5nr identifiers for annotation against Refseq

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash


for FAA in `ls *split* | sed -e 's/ids_md5_subsamples_//g'`
do
echo "#!/bin/bash" > ann_"$FAA".scr
echo "#$ -cwd" >>ann_"$FAA".scr
echo "#$ -j y" >>ann_"$FAA".scr
echo "#$ -S /bin/bash" >>ann_"$FAA".scr

echo python3 m5nr_annotation.py /m5nr/m5nr_refseq.dict ids_md5_subsamples_"$FAA" "$FAA"_out.txt>>ann_"$FAA".scr
chmod +x *.scr; done

#M5nr_refseq.dict is a python dictionary created from the M5nr database
```


```bash
#Concatenate all files from RefSeq
cat split_a* > Refseq_global_ss

#Modify the table 
awk -F "\t" '{print $2"\t"$1"\t"$3}' Refseq_global_ss > Refseq_global_ss2
```


```bash
#Match the output of hitter_table.py with the annotations file

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
perl -e ' $col1=0; $col2=1; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; '  hs_subsamples_global.tsv Refseq_global_ss2 >global_ss_refseq.tsv

```


```bash
#Edit special characters of ontology table
sed -i "s/\//_/g" global_ss_refseq.tsv
sed -i 's/"//g' global_ss_refseq.tsv
sed -i 's/\.//g' global_ss_refseq.tsv
sed -i 's/[//g' global_ss_refseq.tsv
sed -i 's/\[//g' global_ss_refseq.tsv
sed -i 's/\]//g' global_ss_refseq.tsv
sed -i 's/;//g' global_ss_refseq.tsv
sed -i 's/\+//g' global_ss_refseq.tsv
sed -i 's/\#//g' global_ss_refseq.tsv
sed -i 's/\(//g' global_ss_refseq.tsv
sed -i 's/(//g' global_ss_refseq.tsv
sed -i 's/)//g' global_ss_refseq.tsv
sed -i 's/://g' global_ss_refseq.tsv
sed -i 's/&//g' global_ss_refseq.tsv
sed -i 's/*//g' global_ss_refseq.tsv
sed -i 's/=//g' global_ss_refseq.tsv
sed -i 's/?//g' global_ss_refseq.tsv
sed -i 's/^//g' global_ss_refseq.tsv
sed -i 's/\^//g' global_ss_refseq.tsv
sed -i 's/,//g' global_ss_refseq.tsv
sed -i "s/'//g" global_ss_refseq.tsv
sed -i "s/-/_/g" global_ss_refseq.tsv

#Obtain proteins counts table
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33}' global_ss_refseq.tsv > global_ss_refseq_counts.tsv
#Add header to counts table
sed -i '1i id\thydro2\thydro2ss1\thydro2ss2\thydro2ss3\thydro3\thydro3ss1\thydro3ss2\thydro3ss3\thydro4\thydro4ss1\thydro4ss2\thydro4ss3\tsoil1\tsoil1ss1\tsoil1ss2\tsoil1ss3\tsoil2\tsoil2ss1\tsoil2ss2\tsoil2ss3\tsoil3\tsoil3ss1\tsoil3ss2\tsoil3ss3\tsoil4\tsoil4ss1\tsoil4ss2\tsoil4ss3\tsoil5\tsoil5ss1\tsoil5ss2\tsoil5ss3' global_ss_refseq_counts.tsv

#Obtain proteins annotation table.
awk -F "\t" '{print $1"\t"$36}' global_ss_refseq.tsv > global_ss_refseq_onto.tsv
#Add header to annotation table
sed -i '1i id\tprotein1' global_ss_refseq_onto.tsv

```


```bash
#files counts and refseq annotation
global_ss_refseq_onto.tsv
global_ss_refseq_counts.tsv
```

# Add proteins from M5nr that were not found in RefSeq


```bash
#Obtain proteins id in RefSeq
cut -f1 global_ss_refseq_onto.tsv > ref_global_ss_list

#Obtain proteins ids that were not found in RefSeq
cat ids_md5_subsamples.txt ref_global_ss_list | sort | uniq -c | grep '1 '>non_refseq_list 
awk '{print $2}' non_refseq_list>non_refseq_list_globss_ids
```


```bash
#Matching proteins not found in refseq against counts file 

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
for i in `cat non_refseq_list_globss_ids`
do grep -w "$i" hs_subsamples_global.tsv >> globss_norefseq.tsv
done
```


```bash
#Add two columns with increasing numbers of conserved hypothetical proteins as annotations
awk '{print $0"\t""hp_c_"NR-1}' globss_norefseq.tsv> hp_globss_norefseq.tsv 
```

# Add non annotated sequences to analysis


```bash
#Obtain all sequences with hits against M5nr database
for i in `ls *_best.simple.tsv| sed -e 's/_best.simple.tsv//g'`; do awk '{print $1}' "$i"_best.simple.tsv > "$i"_anotados.txt; done 

#Obtain a list of every sequence
for i in `cat muestras`; do grep '>' "$i".faa | sed 's/>//g' | cut -f1 -d ' ' > "$i"_todos.txt; done

#Obtain identifiers of non matching sequences
for i in `cat muestras`; do cat "$i"_anotados.txt "$i"_todos.txt | sort | uniq -c | grep '1 ' | awk '{print $2}' > "$i"_na.txt; done

#Extract non matching sequences.
for i in `cat muestras`; do seqtk subseq "$i".faa "$i"_na.txt > "$i"_na.faa; done

#merge non-matching sequences from all samples
cat *_na.faa > todos_na.faa
```


```bash
#Cluster into gene families at 70% protein identity

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
cd-hit -i todos_na.faa -o todos70 -c 0.70 -n 4 -aL 0.7 -d 0 -M 3000 -T 2 > todos70.cdhit.out

#Convert cdhit output to otu table
perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g;s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' todos70.clstr > todos70.otu
sed -i.bak 's/_/./g;s/RFT./RFT_/g;s/RZPS./RZPS_/g;s/S0./S0_/g' todos70.otu
```


```bash
#Retrieve subsamples proteins 
for HIT in `ls *_ss*.hits | sed -e 's/.hits//g'`
do
awk '{print $2}' "$HIT".hits > "$HIT"_seqname_allhits.txt
done

#Concatenate the id of subsamples proteins with samples unannotated proteins
#Identify duplicates as these are the unannotated subsamples proteins

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
for NA in `ls *na.txt | sed -e 's/_na.txt//g'`
do
cat "$NA"_na.txt "$NA"_ss1_seqname_allhits.txt | sort | uniq -c | grep '2 ' >> "$NA"_ss1.na
cat "$NA"_na.txt "$NA"_ss2_seqname_allhits.txt | sort | uniq -c | grep '2 ' >> "$NA"_ss2.na
cat "$NA"_na.txt "$NA"_ss3_seqname_allhits.txt | sort | uniq -c | grep '2 ' >> "$NA"_ss3.na
done

#Obtain the id unannotated subsamples proteins
for NA in `ls *.na| sed -e 's/.na//g'`
do
awk '{print $2}' "$NA".na > "$NA"_na_seqname.txt
done
```


```bash
#Obtain the file of unannotated proteins from subsamples

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

for i in `ls *.faa | sed -e 's/.faa//g'`
do
seqtk subseq "$i".faa "$i"_ss1_na_seqname.txt > "$i"_ss1_na.faa
seqtk subseq "$i".faa "$i"_ss2_na_seqname.txt > "$i"_ss2_na.faa
seqtk subseq "$i".faa "$i"_ss3_na_seqname.txt > "$i"_ss3_na.faa
done
```


```bash
#Rename subsamples proteins

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
for i in `ls *_na.txt | sed -e 's/_na.txt//g'`
do
sed -i "s/$i/"$i"ss1/g" "$i"_ss1_na.faa
sed -i "s/$i/"$i"ss2/g" "$i"_ss2_na.faa
sed -i "s/$i/"$i"ss3/g" "$i"_ss3_na.faa
done

#Concatenate subsamples proteins
cat *ss*_na.faa > all_ss_na.fas

#Concatenate samples and subsamples proteins
cat *_na.* > ss_todosglobal_na.faa
```


```bash
#Cluster proteins from samples and subsamples into gene families at 70% protein identity

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

cd-hit-2d -i todos70 -i2 all_ss_na.fas -o allss_70 -c 0.70 -n 4 -aL 0.7 -d 0 -M 3000 -T 2 > all_ss_na.cdhit.out
```


```bash
#Convert cdhit output to otu table
perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g;s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' allss_70.clstr > allss_70.otu
sed -i.bak 's/_/./g;s/RFT./RFT_/g;s/RZPS./RZPS_/g;s/S0./S0_/g' allss_70.otu
```


```bash
#rename sequences in subsamples files

for i in `ls *_na.txt | sed -e 's/_na.txt//g'`
do
sed -i "s/$i/"$i"ss1/g" "$i"_ss1.hits
sed -i "s/$i/"$i"ss2/g" "$i"_ss2.hits
sed -i "s/$i/"$i"ss3/g" "$i"_ss3.hits

done

#rename hits files
for i in `ls *_na.txt | sed -e 's/_na.txt//g'`
do
mv "$i"_ss1.hits "$i"ss1.hits
mv "$i"_ss2.hits "$i"ss2.hits
mv "$i"_ss3.hits "$i"ss3.hits
done
```


```bash
#Using hitter_na.py script to obtain abundances of each OTU (hits to M5nr) parsing hits tables and OTU tables. 

#split file in 30 files
split -l 356217 allss_70.otu.bak allss_70.otu.

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

for FAA in `ls *.otu.a* | sed -e 's/allss_70.//g'`
do
echo "#!/bin/bash" > hitna_"$FAA".scr
echo "#$ -cwd" >>hitna_"$FAA".scr
echo "#$ -j y" >>hitna_"$FAA".scr
echo "#$ -S /bin/bash" >>hitna_"$FAA".scr

echo python3 hitter_na.py ss_na_list.txt allss_70."$FAA" ss_na_"$FAA">>hitna_"$FAA".scr
chmod +x *.scr; done
```


```bash
#joining all output files
for i in `ls *.tsv| sed -e 's/.tsv//g'`; do sed '1d' "$i".tsv > "$i".tmp; done
cat *.tmp > hs_na_global_ss.tsv
rm *.tmp

#Add two columns of increasing numbers of hipothetical proteins as annotation
awk '{print $0"\t""hp_"NR-1"\t"}' hs_na_global_ss.tsv >hs_na_global_ss_numerado.tsv
```


```bash
#Unnanotated sequences tables were merged with Refseq annotation and M5nr for complete diversity estimation in the metagenomes

#preparing refseq table
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$36}' global_ss_refseq.tsv > global_ss_refseq_order.tsv

#merge refseq + nonrefseq + hypothetical proteins tables
cat hs_na_global_ss_numerado.tsv global_ss_refseq_order.tsv hp_globss_norefseq.tsv  > all_proteins_global_ss.tsv

#proteins counts table
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33}' all_proteins_global_ss.tsv > all_proteins_global_ss_counts.tsv
#proteins annotation table.
awk -F "\t" '{print $1"\t"$34}' all_proteins_global_ss.tsv > all_proteins_global_ss_onto.tsv

```


```bash
#add header to counts table
sed -i '1i Identifier\thydro2\thydro2ss1\thydro2ss2\thydro2ss3\thydro3\thydro3ss1\thydro3ss2\thydro3ss3\thydro4\thydro4ss1\thydro4ss2\thydro4ss3\tsoil1\tsoil1ss1\tsoil1ss2\tsoil1ss3\tsoil2\tsoil2ss1\tsoil2ss2\tsoil2ss3\tsoil3\tsoil3ss1\tsoil3ss2\tsoil3ss3\tsoil4\tsoil4ss1\tsoil4ss2\tsoil4ss3\tsoil5\tsoil5ss1\tsoil5ss2\tsoil5ss3' all_proteins_global_ss_counts.tsv
#add header to annotation table
sed -i '1i Identifier\tprotein1' all_proteins_global_ss_onto.tsv
#create another colum for onthology
awk '{print $1"\t"$2"\t"$2}' all_proteins_global_ss_onto.tsv > all_proteins_global_ss_onto2.tsv
#Edit special characters of ontology table
sed -i "s/\//_/g" all_proteins_global_ss_onto2.tsv
sed -i 's/"//g' all_proteins_global_ss_onto2.tsv
sed -i 's/\.//g' all_proteins_global_ss_onto2.tsv
sed -i 's/[//g' all_proteins_global_ss_onto2.tsv
sed -i 's/\[//g' all_proteins_global_ss_onto2.tsv
sed -i 's/\]//g' all_proteins_global_ss_onto2.tsv
sed -i 's/;//g' all_proteins_global_ss_onto2.tsv
sed -i 's/\+//g' all_proteins_global_ss_onto2.tsv
sed -i 's/\#//g' all_proteins_global_ss_onto2.tsv
sed -i 's/\(//g' all_proteins_global_ss_onto2.tsv
sed -i 's/(//g' all_proteins_global_ss_onto2.tsv
sed -i 's/)//g' all_proteins_global_ss_onto2.tsv
sed -i 's/://g' all_proteins_global_ss_onto2.tsv
sed -i 's/&//g' all_proteins_global_ss_onto2.tsv
sed -i 's/*//g' all_proteins_global_ss_onto2.tsv
sed -i 's/=//g' all_proteins_global_ss_onto2.tsv
sed -i 's/?//g' all_proteins_global_ss_onto2.tsv
sed -i 's/^//g' all_proteins_global_ss_onto2.tsv
sed -i 's/\^//g' all_proteins_global_ss_onto2.tsv
sed -i 's/,//g' all_proteins_global_ss_onto2.tsv
sed -i "s/'//g" all_proteins_global_ss_onto2.tsv
sed -i "s/-/_/g" all_proteins_global_ss_onto2.tsv
all_proteins_global_ss_onto2.tsv >all_proteins_global_ss_onto.tsv
```


```bash
#Final files
all_proteins_global_ss_counts.tsv
all_proteins_global_ss_onto.tsv
```

# Pangenome construction and metagenome search


```bash
#Creating gff3 files with Prokka
#This was done for every genome of every bacterial genus

prokka --kingdom Bacteria --outdir prokka_Strain1 --genus Luteolibacter --locustag Strain1 $PATH/Strain1_genomic.fna

#Running Roary

bin/roary -f output -e -n -i 50 -v ./*.gff
anaconda3/bin/FastTreeMP -nt -gtr core_gene_alignment.aln > tree_core.nwk
python roary_plots.py ./output_outgroup/tree_core.nwk ./output_outgroup/gene_presence_absence.csv
Rscript create_pan_genome_plots.R

#Running promer for every pangenome reference using the concatenated reads of the metagenomes as query (allhydro_hb100_contigs.fasta.numbered.fas)#

promer pan_genome_reference.fa allhydro_hb100_contigs.fasta.numbered.fas 
mv out.delta luteo.delta; done

#To edit promer graphs

show-coords -c -d -k -l -r -T luteo.delta > luteo.show-coords.tsv
python3 promer_deid_v9gmv.py luteo.show-coords.tsv pangenome_Luteolibacter;
```

# Taxonomic assignment of contigs


```bash
#crear un archivo con la lista de los nombres 

#!/bin/bash
#bash
SEQS=hybrid_ass/assh_merge
DB=dbs/kraken2
BIN=binc/kraken2
OUT=protein_tax/kraken

COUNT=0

for FAA in `cat list`
do let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo "$BIN"/kraken2 --use-names --report "$OUT/$FAA".tax --use-mpa-style --db "$DB" "$SEQS/$FAA"_hb100_contigs.fasta.numbered.fas ">" "$OUT/$FAA".kout >> $*.$COUNT.scr
chmod +x *.scr; done

#list file contains the names of the samples
```


```bash
#Concatenate all the results into a single file
cat *.kout > all_hydro_soil.kout
#cleaning table
awk -F "\t" '{print $2"\t"$3}' all_hydro_soil.kout > hs_krkn.tmp
sed -i "s/taxid //g" hs_krkn.tmp
sed -i 's/ (/\t/g' hs_krkn.tmp
sed -i 's/)//g' hs_krkn.tmp
sed -i 's/ /_/g' hs_krkn.tmp

```

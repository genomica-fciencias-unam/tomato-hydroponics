```bash

# Experimenting with PANDASEQ
# Manual trimming based on quality is required. This protocol was effective with MiSEQ using 2x250 reads.

# Assemble with PANDASEQ
ls *.fastq.gz | perl -pe 's/_R.*.fastq.gz//g' | sort | uniq >list # Generate the list
ln -s scripts/assemblyPANDA.sh assemblyPANDA.sh # Create a symbolic link to the PANDASEQ assembly script
bash assemblyPANDA.sh JOB_NAME # Run the script that generates assembly jobs using PANDASEQ (pre-requisite)
for N in `ls *.scr`; do qsub $N; done # Submit the assembly jobs to the cluster

# Rename sequences with simple identifiers (see label section)
# Rename sequences with identifiers using the sequence name and number them:
# >Sequence_0 
# ATTAC
# >Sequence_1
# TTCAT
# >Sequence_2
# CCTAT

# Perform this step for each sample with individual labels per sample
perl scripts/header.fasta.numbers.pl PREFIX file_name.fasta 

# Concatenate all samples from the study
cat *.numbered.fasta >complete_study.fas

# Group sequences at 97% identity
# If using the deep-thought cluster, use the parallelization script
# If using the cuallicua cluster, use the script directly but submit via the queue manager (qsub)

# On deep-thought:
ln -s scripts/cd-clust.sh . 
bash cd-clust.sh complete_study.fas

# On cuallicua:
qsub -N JOB_NAME -b y -j y -cwd -V "cd-hit-est -c 0.97 -T 25 -M 0 -i complete_study.fas -o output.clstr"

# Generate OTU table
perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g; s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"' ; exit}' output.clstr.clstr >mp_16s.otu



# Command to remove singletons from the OTU file
awk 'NF < 3' mp_16S.otu | cut -f1 >ids_REMOVE_singletons.txt
wc -l ids_REMOVE_singletons.txt 
cut -f1 mp_16S.otu >rawids

# Generate a file with the identities of singletons
awk 'NR==FNR { values[$0]; next } !($0 in values)' ids_REMOVE_singletons.txt rawids  >repsetWOsingletons.txt

# Extract FASTA without singletons 
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV'repsetWOsingletons.txt mp_16S.rep.fna >lda_nosingletons.fna
grep -c ">" lda_nosingletons.fna

# First search for contaminants (non-16S sequences)
parallel_assign_taxonomy_blast.py -i lda_nosingletons.fna -o ldaNO16S_screen1 -r /qiime/gg_otus-13_8-release/rep_set/70_otus.fasta -t /qiime/gg_otus-13_8-release/taxonomy/70_otu_taxonomy.txt -X LDA1

# List sequences with no hits against the Greengenes 70 database
cat ldaNO16S_screen1/lda_nosingletons_tax_assignments.txt | grep "No blast hit" | cut -f1 >idsnohit.txt
wc -l idsnohit.txt

# Extract FASTA of sequences with no hit
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' idsnohit.txt lda_nosingletons.fna  >nohitGG.fasta

# Verify sequences with no hits in the latest Silva database
parallel_assign_taxonomy_blast.py -i nohitGG.fasta -o LDAsinhitGG -r /databases/silva/qiime/dada_silva138.fna -t /databases/silva/qiime/dada_silva138.tax.txt

# Obtain identifiers for sequences without singletons and contaminants (screened 16S)
grep -v "No blast hit" ldaNO16S_screen1/lda_nosingletons_tax_assignments.txt | cut -f1 >ggret1.txt
grep -v "No blast hit" LDAsinhitGG/nohitGG_tax_assignments.txt | cut -f1 >ggret2.txt
cat ggret*.txt >ret_screened_nosingleton.txt

wc -l ret_screened_nosingleton.txt 

# Extract FASTA of sequences with no hit
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ret_screened_nosingleton.txt  lda_nosingletons.fna >nosingleton_screened.fna

# Align to check for chimeras:
parallel_align_seqs_pynast.py -i nosingleton_screened.fna -o chimera_align -O 400 -X chimx

# Search for chimeras:
parallel_identify_chimeric_seqs.py -m blast_fragments -i nosingleton_screened.fna -a chimera_align/nosingleton_screened_aligned.fasta -o singleton_chimera.txt -X chimerablast --id_to_taxonomy_fp /qiime/gg_otus-13_8-release/taxonomy/97_otu_taxonomy.txt -r /qiime/gg_otus-13_8-release/rep_set/97_otus.fasta -O 300
wc -l singleton_chimera.txt

# Obtain chimera IDs
awk '{print $1}' singleton_chimera.txt >chimeraids

# Remove these identities from the second file using the first file
awk 'NR==FNR { values[$0]; next } !($0 in values)' chimeraids ret_screened_nosingleton.txt >ids_nochimera_screened_nosingleton.txt

wc -l ids_nochimera_screened_nosingleton.txt 

# Generate a new OTU table without chimeras, filtered, without singletons:
awk 'NR==FNR { values[$0]; next } $1 in values { print }' ids_nochimera_screened_nosingleton.txt mp_16S.otu >nochimera_screened_singleton.otu

# New FASTA without chimeras, filtered, without singletons:
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids_nochimera_screened_nosingleton.txt nosingleton_screened.fna >nochimera_screened_nosingleton.fna
grep -c ">" nochimera_screened_nosingleton.fna 

# Assign taxonomy
parallel_assign_taxonomy_blast.py -i nochimera_screened_nosingleton.fna -o tax -r /databases/silva/qiime/dada_silva138.fna -t /databases/silva/qiime/dada_silva138.tax.txt -O 300 -X tom

# Clean taxonomy
awk '{ for(i=1; i<=NF-2; i++) printf "%s%s", $i, (i<NF-2) ? OFS : ORS }' tax/nochimera_screened_nosingleton_tax_assignments.txt | perl -pe 's/\ /\_/g;s/_k/\tk/;s/\;_/\t/g' >tax_clean.txt

# Remove mitochondria and chloroplasts:
grep -i 'chloroplast' tax_clean.txt | cut -f1 >tmprm1
grep -i 'mitochondria' tax_clean.txt | cut -f1 >>tmprm1

# Generate a new OTU table without chimeras, filtered, without singletons, without mitochondria and chloroplasts:
awk 'NR==FNR { values[$0]; next } !($0 in values)' tmprm1 ids_nochimera_screened_nosingleton.txt >nochimera_screened_singleton_cpmt.ids
awk 'NR==FNR { values[$0]; next } $1 in values { print }' nochimera_screened_singleton_cpmt.ids nochimera_screened_singleton.otu >nochimera_screened_singleton_cpmt.otu
wc -l nochimera_screened_singleton_cpmt.otu 

# FASTA of OTUs without chimeras, filtered, without singletons, without mitochondria and chloroplasts:
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV'  nochimera_screened_singleton_cpmt.ids nochimera_screened_nosingleton.fna >nochimera_screened_singleton_cpmt.fna

# Generate alignment of filtered OTUs:
ssu-prep -f nochimera_screened_singleton_cpmt.fna tomatotax 40 prefixcluster.txt
bash tomatotax.ssu-align.sh
# Once it finishes in the cluster
ssu-mask --stk2afa tomatotax

# Clean alignment
perl -i.bak -pe 's/\./-/g;s/U/T/gi' tomatotax/tomatotax.bacteria.afa

# Build Phylogeny
./fasttreeMP -nt -fastest ssu-align/tomatotax/tomatotax.bacteria.afa >tomato.nwk

#Another way to create the OTU table:
#Create a list of samples: 
grep ">" mp_16S.fasout | perl -pe 's/\>//g;s/\_/\t/g' | cut -f1 | sort | uniq >samples
python3 otu_tab.py samples nochimera_screened_singleton_cpmt.otu

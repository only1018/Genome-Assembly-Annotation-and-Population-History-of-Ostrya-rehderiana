#Primary genome assembly
hifiasm -o 101-Ore-le-1 -t 30 `ls ./*.fasta.gz` \ 
        --h1 FDHC250018947-1b_L1_1.clean.rd.fq.gz FDHC250018947-1b_L2_1.clean.rd.fq.gz FDHC250018947-1b_L3_1.clean.rd.fq.gz \
        --h1 FDHC250018947-1b_L1_2.clean.rd.fq.gz FDHC250018947-1b_L2_2.clean.rd.fq.gz FDHC250018947-1b_L3_2.clean.rd.fq.gz

#Read mapping-based assembly validation
bwa index -a bwtsw genome.fa
bwa mem -k 32 -w 10 -B 3 -O 11 -E 4 -t 8 genome.fa *.clean.rd.fq.gz | samtools view -b -S - > *.clean.rd.fq.bam

#Genome completeness assessment
busco -i genome.fa -l embryophyta_odb10 -c 15

#Consensus quality (QV) estimation
merqury  19 new_genome.fa input.fofn


#Hi-C assisted genome scaffolding
bowtie2-build --threads 8 genome.fa ref
HiC-Pro -v -S -t 100 -m 100000 -s 100 -l 800 -a -d 1000 -f dpnii_contigs.bed -r all.bam -o hicpro
juicer -a genome.assembly -c 8 -i merged_nodups.txt -o ./

#chloroplast and mitochondrial genome assembly
oatk -k 1001 -c 150 -t 8 -m angiosperm_mito.fam -o oatk $mark_hifi.fasta.gz

#Filtering of organellar sequences
blastall -p blastn -i ./*.fasta -d *.fasta -o blast.m8 -m 8 -a 8 -e 1e-5

#Gap filling
tgsgapcloser --racon \
             --tgstype ont \
             --minmap_arg ' -x ava-ont' \
             --min_nread 5

#Telomere identification
tidk search --string AAACCCT --output test2 --dir $mark --window 100000 $genome

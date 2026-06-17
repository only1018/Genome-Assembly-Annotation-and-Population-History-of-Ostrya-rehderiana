#Primary genome assembly
hifiasm -o 101-Ore-le-1 -t 30 `ls ./*.fasta.gz` \ 
        --h1 FDHC250018947-1b_L1_1.clean.rd.fq.gz FDHC250018947-1b_L2_1.clean.rd.fq.gz FDHC250018947-1b_L3_1.clean.rd.fq.gz \
        --h1 FDHC250018947-1b_L1_2.clean.rd.fq.gz FDHC250018947-1b_L2_2.clean.rd.fq.gz FDHC250018947-1b_L3_2.clean.rd.fq.gz

#Read mapping-based assembly validation
Bwa_pipeline.pl genome.fa ./BWA reads.list

#Genome completeness assessment
busco -i genome.fa -l embryophyta_odb10 -c 15

#Consensus quality (QV) estimation
merqury  19 new_genome.fa input.fofn


#Hi-C assisted genome scaffolding
hic.pl \
  -fasta genome.fa \
  -reads reads.list \
  -chr_num 8 \
  -sample 101-Ore-le-1 \
  -class plant \
  -enz DpnII \
  -workdir ./ \
  -sge_queue denovo.q

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

#Read mapping
bwa mem -t 20 reference_genome.fa reads_R1.fq.gz reads_R2.fq.gz | \
samtools view -bS - > sample.bam

#Sorting and indexing
samtools sort -@ 20 -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam

#Variant calling and diploid consensus generation
samtools mpileup -C50 -uf reference_genome.fa sample.sorted.bam | \
bcftools call -c - | \
vcfutils.pl vcf2fq -d 25 -D 152 > sample.psmc.fq

#Convert consensus to PSMC format
fq2psmcfa -q20 sample.psmc.fq > sample.psmcfa

#Run PSMC
psmc -N30 -t5 -r5 -p "1+1+1+1+25*2+4+6" -o sample.psmc sample.psmcfa

#Bootstrap replicates (robustness test)
splitfa sample.psmcfa > split.psmcfa

seq 100 | xargs -i echo "psmc -N30 -t5 -r5 -b -p \"1+1+1+1+25*2+4+6\" \
-o bootstrap_{}.psmc split.psmcfa" | sh

#Plot results
psmc_plot.pl -u 2.182e-8 -g 10 sample_plot sample.psmc

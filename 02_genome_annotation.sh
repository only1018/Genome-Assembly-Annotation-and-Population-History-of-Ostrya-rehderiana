#Repeat Annotation
RepeatModeler -database model_genome -engine ncbi -threads 30 -LTRStruct
RepeatMasker -a -no_is -norna -pa 30 -lib genome.libraries.fa -q genome.fa -gff

#Non-coding RNA Annotation
blastall -p blastn -e 1e-10 -v 10000 -b 10000 -d genome.fa -i plant_rRNA.fa -o genome.fa.rRNA.blast
cmsearch --cpu 20 db.1 genome.fa > genome.fa.snRNA.db.cmsearch
tRNAscan-SE -Q genome.fa -o genome.fa.tRNA -f genome.fa.tRNA.structure

#Transcriptome-based Annotation
Trinity --seqType fq \
  --left *_1.clean.fq.gz \
  --right *_2.clean.fq.gz \
  --CPU 32 \
  --max_memory 75G \
  --normalize_reads \
  --full_cleanup \
  --min_glue 2 \
  --min_kmer_cov 2 \
  --output trinity.out

hisat2 -p 16 -x hisat_genome.fa -1 *_1.clean.fq.gz -2 *_2.clean.fq.gz
stringtie *.sort.bam -p 16 -o out.gtf

#De novo Annotation
augustus --AUGUSTUS_CONFIG_PATH=conf \
    --species=pasa1 --uniqueGeneId=true \
    --noInFrameStop=true --gff3=on --genemodel=complete \
    --strand=both genome_masked.fa > genome_masked.fa.Augustus

snap genome.fa -gff genome.train.nr.uniq.filter.gff

#Homology-based Annotation
tblastn -query homo.pep -db genome.chunk.fa -outfmt 6 -evalue 1e-05
homo_solar -a prot2genome2 -z -f m8 *.pep.blast > *.pep.blast.solar
genewise --gw *.gw --pep *.filter.pep --type 1 >out.aln
genewise --gw *.gw --pep *.filter.pep --type 2 >out.mut

#Evidence Integration of Gene Annotations
evm_partition_inputs --outputdir $PWD \
    --genome $PWD/*.fa \
    --gene_predictions $PWD/denovo.gff \
    --protein_alignments $PWD/homolog.gff \
    --transcript_alignments $PWD/transcripts.gff \
    --segmentSize 200000 \
    --overlapSize 20000 \
    --min_intron_length 20 \
    --partition_listing evm.partition

#Functional Annotation
interproscan -cpu 20 -format tsv -iprlookup -dp \
    -goterms -i *.fa -o out.fa.iprscan
blastall -b 500 -v 500 -F F -p blastp -e 1e-4 \
-i pasa2.longest.fa -o pasa2.longest.pep.fa.uniprot_sprot.fasta.blast \
-d swissprot NR KEGG

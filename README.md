# Genome Assembly, Annotation and Population History of *Ostrya rehderiana*
This pipeline outlines the generation of a genome assembly, annotation, and inference of effective population size ($N_\mathrm{e}$) dynamics for the critically endangered species *Ostrya rehderiana*.
- 01_genome_assembly.sh outlines the genome assembly workflow;
- 02_genome_annotation.sh outlines the genome annotation pipeline；
- 03_centromere_prediction.sh provides the pipeline to predict putative centromeric regions, and validate their assembly continuity using ONT long reads.
- 04_psmc.sh provides the pipeline to get inference of effective population size dynamics via PSMC；
- 05_psmc_plot.R provides scripts for visualization of PSMC-derived demographic trajectories.

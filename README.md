# Analysis of Hawaiian Coral Genomes

**Description of directories**

- `04_Stats_Original_Montipora_capitata_Assembly/` Scripts and results for BUSCO analysis of original *M. capitata* assembly (Shumaker, A., Putnam, H.M., Qiu, H. et al. Genome analysis of the rice coral Montipora capitata. Sci Rep 9, 2571 (2019). https://doi.org/10.1038/s41598-019-39274-3)
- `01_ExtraAnalysis-Montipora_capitata_KBHIv3/01_explore_extrachromosomal_scaffolds/` Scripts and results for comparison between chromosomal and non-chromosomal contigs from the new *M. capitata* assembly. The `blastn_across_scaffold_sets/` sub-directory contains a README.md describing the comparison between the contigs using BLASTN and the `blastp_across_protein_sets` sub-directory contains a README.md describing the comparison between the proteins redicted on the two different sets of contigs using BLASTP.
- `01_ExtraAnalysis-Montipora_capitata_KBHIv3/04_BUSCO_ScaffoldsSetsSplit/` Scripts and results for BUSCO analysis of the new *M. capitata* assembly but split into chromosomal (the `genome.fa` and `pep.faa` files) and non-chromosomal sets (the `genome_extra.fa` and `pep_extra.faa` files). This directory also contains a R markdown file (`BUSCO_results.Rmd`) which will take the results from BUSCO nd will produce a combined plot showing the different BUSCO sets for each of the analyzed genomes (both generated by us and published).
- `01_ExtraAnalysis-Montipora_capitata_KBHIv3/05_GeneStats_ScaffoldsSetsSplit` R markdown script (`Stats2table.Rmd`) for generating an excel table with all of the basic genome assembly and predicted gene stats (e.g., N50, no. contigs, No. proteins, etc.) for each of the analyzed genomes (both generated by us and published).
- `shared_scripts/` This directory contains copies of generic scripts that I use for text parsing or processing of fasta or table files. NOTE: In the \*.md files in the other directories these scripts are generally shown to be in a higher level directory called `02_Scripts`.


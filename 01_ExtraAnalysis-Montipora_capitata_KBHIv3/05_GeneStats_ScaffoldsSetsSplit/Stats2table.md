---
title: "Plot QC results for each genome"
author: "Timothy Stephens"
date: "23/02/2022"
output: 
  html_document:
    keep_md: yes
---



# Setup

Setup R env. Load packages and set default image export formats, size and resolution.


```r
knitr::opts_chunk$set(echo = TRUE,
                      fig.height = 8, 
                      fig.width = 12, 
                      dev = c("png", "pdf"),
                      dpi = 1000)
library(readxl)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(tibble)
library(ggplot2)
library(reshape2)
library(scales)
library(ggdendro)
library(cowplot)
library(RColorBrewer)
library(phylogram)
library(gridExtra)
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```

```r
library(ggrepel)
library(writexl)
library(stringr)
options(scipen = 999) #Prevent scientific notation
```


## Color palets

Color palet for BUSCO gene categories.


```r
my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442")
```





# Format Gene Stats info


```bash
i=1
for F in ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/Montipora_*/01_stats/*.GeneStats.tsv \
         ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/Pocillopora_*/01_stats/*.GeneStats.tsv \
         Montipora_capitata_KBHIv3.chromosomes.GeneStats.tsv \
         Montipora_capitata_KBHIv3.extrachromosomal.GeneStats.tsv;
do 
  filelist[i]="$F"
  i=$i+1
done

awk -F'\t' 'BEGIN {L="SampleID"} {L=L"\t"$1} END {print L} ' "${filelist[1]}" > all_genomes-01_stats-results.tsv
for F in "${filelist[@]}";
do
  awk -F'\t' 'BEGIN {L=""} { if(NR==1) {X=$2; gsub(".assembly.fasta", "", $2); L=L"\t"$2"\t"X} else {L=L"\t"$2}} END {print L} ' "${F}" | cut -f2- >> all_genomes-01_stats-results.tsv
done
```





# Format BUSCO info

Combine BUSCO gene count and percent into a single cell and format for inport into R.


```bash
i=1
for F in ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/Montipora_*/02_busco \
         ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/Pocillopora_*/02_busco;
do 
  filelist[i]="$F"
  i=$i+1
done

## Genome - eukaryota_odb10
awk -F'\t' 'BEGIN {L="SampleID"} {L=L"\t"$1} END {print L} ' "${filelist[1]}/genome.fa.busco_eukaryota_odb10.results.txt" \
  > all_genomes-02_busco-Genome_eukaryota_odb10-results4table.tsv
for F in "${filelist[@]}";
do
  P=$(echo "${F}" | awk -F'/' '{print $10}')
  awk -F'\t' -v P="$P" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' "${F}/genome.fa.busco_eukaryota_odb10.results.txt" \
    >> all_genomes-02_busco-Genome_eukaryota_odb10-results4table.tsv
done
awk -F'\t' -v P="Montipora_capitata_KBHIv3.chromosomes" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' \
  ../04_BUSCO_ScaffoldsSetsSplit/genome.fa.busco_eukaryota_odb10.results.txt \
    >> all_genomes-02_busco-Genome_eukaryota_odb10-results4table.tsv
awk -F'\t' -v P="Montipora_capitata_KBHIv3.extrachromosomal" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' \
  ../04_BUSCO_ScaffoldsSetsSplit/genome_extra.fa.busco_eukaryota_odb10.results.txt \
    >> all_genomes-02_busco-Genome_eukaryota_odb10-results4table.tsv


## Genome - metazoa_odb10
awk -F'\t' 'BEGIN {L="SampleID"} {L=L"\t"$1} END {print L} ' "${filelist[1]}/genome.fa.busco_metazoa_odb10.results.txt" \
  > all_genomes-02_busco-Genome_metazoa_odb10-results4table.tsv
for F in "${filelist[@]}";
do
  P=$(echo "${F}" | awk -F'/' '{print $10}')
  awk -F'\t' -v P="$P" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' "${F}/genome.fa.busco_metazoa_odb10.results.txt" \
    >> all_genomes-02_busco-Genome_metazoa_odb10-results4table.tsv
done
awk -F'\t' -v P="Montipora_capitata_KBHIv3.chromosomes" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' \
  ../04_BUSCO_ScaffoldsSetsSplit/genome.fa.busco_metazoa_odb10.results.txt \
    >> all_genomes-02_busco-Genome_metazoa_odb10-results4table.tsv
awk -F'\t' -v P="Montipora_capitata_KBHIv3.extrachromosomal" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' \
  ../04_BUSCO_ScaffoldsSetsSplit/genome_extra.fa.busco_metazoa_odb10.results.txt \
    >> all_genomes-02_busco-Genome_metazoa_odb10-results4table.tsv

## Protein - eukaryota_odb10
awk -F'\t' 'BEGIN {L="SampleID"} {L=L"\t"$1} END {print L} ' "${filelist[1]}/pep.faa.busco_eukaryota_odb10.results.txt" \
  > all_genomes-02_busco-Protein_eukaryota_odb10-results4table.tsv
for F in "${filelist[@]}";
do
  P=$(echo "${F}" | awk -F'/' '{print $10}')
  awk -F'\t' -v P="$P" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' "${F}/pep.faa.busco_eukaryota_odb10.results.txt" \
    >> all_genomes-02_busco-Protein_eukaryota_odb10-results4table.tsv
done
awk -F'\t' -v P="Montipora_capitata_KBHIv3.chromosomes" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' \
  ../04_BUSCO_ScaffoldsSetsSplit/pep.faa.busco_eukaryota_odb10.results.txt \
    >> all_genomes-02_busco-Protein_eukaryota_odb10-results4table.tsv
awk -F'\t' -v P="Montipora_capitata_KBHIv3.extrachromosomal" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' \
  ../04_BUSCO_ScaffoldsSetsSplit/pep_extra.faa.busco_eukaryota_odb10.results.txt \
    >> all_genomes-02_busco-Protein_eukaryota_odb10-results4table.tsv

## Protein - metazoa_odb10
awk -F'\t' 'BEGIN {L="SampleID"} {L=L"\t"$1} END {print L} ' "${filelist[1]}/pep.faa.busco_metazoa_odb10.results.txt" \
  > all_genomes-02_busco-Protein_metazoa_odb10-results4table.tsv
for F in "${filelist[@]}";
do
  P=$(echo "${F}" | awk -F'/' '{print $10}')
  awk -F'\t' -v P="$P" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' "${F}/pep.faa.busco_metazoa_odb10.results.txt" \
    >> all_genomes-02_busco-Protein_metazoa_odb10-results4table.tsv
done
awk -F'\t' -v P="Montipora_capitata_KBHIv3.chromosomes" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' \
  ../04_BUSCO_ScaffoldsSetsSplit/pep.faa.busco_metazoa_odb10.results.txt \
    >> all_genomes-02_busco-Protein_metazoa_odb10-results4table.tsv
awk -F'\t' -v P="Montipora_capitata_KBHIv3.extrachromosomal" 'BEGIN {L=P} {L=L"\t"$3" ("$2")"} END {print L} ' \
  ../04_BUSCO_ScaffoldsSetsSplit/pep_extra.faa.busco_metazoa_odb10.results.txt \
    >> all_genomes-02_busco-Protein_metazoa_odb10-results4table.tsv
```





# Load data


```r
genome.stats <- read.table("all_genomes-01_stats-results.tsv", sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
genome.eukaryota <- read.table("all_genomes-02_busco-Genome_eukaryota_odb10-results4table.tsv", sep='\t', header=TRUE, check.names=FALSE)
genome.metazoa <- read.table("all_genomes-02_busco-Genome_metazoa_odb10-results4table.tsv", sep='\t', header=TRUE, check.names=FALSE)
protein.eukaryota <- read.table("all_genomes-02_busco-Protein_eukaryota_odb10-results4table.tsv", sep='\t', header=TRUE, check.names=FALSE)
protein.metazoa <- read.table("all_genomes-02_busco-Protein_metazoa_odb10-results4table.tsv", sep='\t', header=TRUE, check.names=FALSE)
```





# Make combined table

Extract just the columns (in a logical order) from each datasets to output into an excel sheet for publication. 
Add extra empty columns which will separate the data from different sections (saves having to add this formatting later by hand). 
Also add "XX--" prefix to BUSCO columns names to make then unique when joining the tables. Will remove later.


```r
t1 <- genome.stats %>%
  filter(SampleID!="Pocillopora_meandrina_KBHIv1") %>%
  select("SampleID", 
         "Total scaffold length (bp)", 
         "Total contig length (bp)",
         "Number of scaffolds",
         "Number of contigs",
         "N50 of scaffolds (bp)",
         "N50 of contigs (bp)",
         "Percent gaps") %>%
  mutate("Genome assembly"="") %>%
  relocate(last_col(), .before = "Total scaffold length (bp)") %>%
  mutate(across(c("Total scaffold length (bp)", 
                  "Total contig length (bp)",
                  "Number of scaffolds",
                  "Number of contigs",
                  "N50 of scaffolds (bp)",
                  "N50 of contigs (bp)"), comma))

t2 <- genome.metazoa %>%
  filter(SampleID!="Pocillopora_meandrina_KBHIv1") %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"genome\" completeness (v5.0; metazoa_odb10; 954 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("GM--Complete (no.)" = "Complete") %>%
  rename("GM--  Single-copy (no.)" = "Single-copy") %>%
  rename("GM--  Duplicated (no.)" = "Duplicated") %>%
  rename("GM--Fragmented (no.)" = "Fragmented") %>%
  rename("GM--Missing (no.)" = "Missing")
  
t3 <- genome.eukaryota %>%
  filter(SampleID!="Pocillopora_meandrina_KBHIv1") %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"genome\" completeness (v5.0; eukaryota_odb10; 255 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("GE--Complete (no.)" = "Complete") %>%
  rename("GE--  Single-copy (no.)" = "Single-copy") %>%
  rename("GE--  Duplicated (no.)" = "Duplicated") %>%
  rename("GE--Fragmented (no.)" = "Fragmented") %>%
  rename("GE--Missing (no.)" = "Missing")
 
t4 <- genome.stats %>%
  filter(SampleID!="Pocillopora_meandrina_KBHIv1") %>%
  select("SampleID", 
         "Number of protein-coding genes", 
         "Average protein-coding gene length (CDS+introns; bp)",
         "Average transcript (CDS) length (bp)",
         "Average number of CDS per gene",
         "Average CDS length (bp)",
         "Number of single-CDS genes",
         "Percent single-CDS genes",
         "Average intron length (between CDS; bp)") %>%
  mutate("Protein-coding genes"="") %>%
  relocate(last_col(), .before = "Number of protein-coding genes") %>%
  mutate("Protein-coding transcripts"="") %>%
  relocate(last_col(), .before = "Average transcript (CDS) length (bp)") %>%
  mutate("Single-CDS genes"="") %>%
  relocate(last_col(), .before = "Number of single-CDS genes") %>%
  mutate("Introns between CDS"="") %>%
  relocate(last_col(), .before = "Average intron length (between CDS; bp)") %>%
  mutate(across(c("Number of protein-coding genes", 
                  "Average protein-coding gene length (CDS+introns; bp)",
                  "Average transcript (CDS) length (bp)",
                  "Average number of CDS per gene",
                  "Average CDS length (bp)",
                  "Number of single-CDS genes",
                  "Average intron length (between CDS; bp)"), comma))

t5 <- protein.metazoa %>%
  filter(SampleID!="Pocillopora_meandrina_KBHIv1") %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"protein\" completeness (v5.0; metazoa_odb10; 954 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("PM--Complete (no.)" = "Complete") %>%
  rename("PM--  Single-copy (no.)" = "Single-copy") %>%
  rename("PM--  Duplicated (no.)" = "Duplicated") %>%
  rename("PM--Fragmented (no.)" = "Fragmented") %>%
  rename("PM--Missing (no.)" = "Missing")
  
t6 <- protein.eukaryota %>%
  filter(SampleID!="Pocillopora_meandrina_KBHIv1") %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"protein\" completeness (v5.0; eukaryota_odb10; 255 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("PE--Complete (no.)" = "Complete") %>%
  rename("PE--  Single-copy (no.)" = "Single-copy") %>%
  rename("PE--  Duplicated (no.)" = "Duplicated") %>%
  rename("PE--Fragmented (no.)" = "Fragmented") %>%
  rename("PE--Missing (no.)" = "Missing")

t <- merge(t1, t2, by="SampleID", sort=FALSE, no.dups=FALSE)
t <- merge(t,  t3, by="SampleID", sort=FALSE, no.dups=FALSE)
t <- merge(t,  t4, by="SampleID", sort=FALSE, no.dups=FALSE)
t <- merge(t,  t5, by="SampleID", sort=FALSE, no.dups=FALSE)
t <- merge(t,  t6, by="SampleID", sort=FALSE, no.dups=FALSE)
t <- t %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(rowname = str_remove(rowname, "GM--")) %>%
  mutate(rowname = str_remove(rowname, "GE--")) %>%
  mutate(rowname = str_remove(rowname, "PM--")) %>%
  mutate(rowname = str_remove(rowname, "PE--"))
```



```r
write_xlsx(t, "all_genomes-Combined-results.xlsx")
```





# Session Info


```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-apple-darwin18.6.0 (64-bit)
## Running under: macOS Mojave 10.14.6
## 
## Matrix products: default
## BLAS/LAPACK: /usr/local/Cellar/openblas/0.3.7/lib/libopenblasp-r0.3.7.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] stringr_1.4.0      writexl_1.4.0      ggrepel_0.9.1      gridExtra_2.3     
##  [5] phylogram_2.1.0    RColorBrewer_1.1-2 cowplot_1.1.1      ggdendro_0.1.23   
##  [9] scales_1.1.1       reshape2_1.4.4     ggplot2_3.3.5      tibble_3.1.6      
## [13] dplyr_1.0.8        readxl_1.3.1      
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.1.2 xfun_0.30        bslib_0.3.1      purrr_0.3.4     
##  [5] lattice_0.20-45  colorspace_2.0-3 vctrs_0.3.8      generics_0.1.2  
##  [9] htmltools_0.5.2  yaml_2.3.5       utf8_1.2.2       rlang_1.0.2     
## [13] jquerylib_0.1.4  pillar_1.7.0     glue_1.6.2       withr_2.5.0     
## [17] DBI_1.1.2        lifecycle_1.0.1  plyr_1.8.6       munsell_0.5.0   
## [21] gtable_0.3.0     cellranger_1.1.0 evaluate_0.15    knitr_1.37      
## [25] fastmap_1.1.0    parallel_3.6.1   fansi_1.0.2      Rcpp_1.0.8      
## [29] jsonlite_1.8.0   digest_0.6.29    stringi_1.7.6    grid_3.6.1      
## [33] cli_3.2.0        tools_3.6.1      magrittr_2.0.2   sass_0.4.0      
## [37] crayon_1.5.0     ape_5.6-1        pkgconfig_2.0.3  ellipsis_0.3.2  
## [41] MASS_7.3-55      assertthat_0.2.1 rmarkdown_2.12   rstudioapi_0.13 
## [45] R6_2.5.1         nlme_3.1-155     compiler_3.6.1
```

---
title: "Plot BUSCO results"
author: "Timothy Stephens"
date: "07/06/2022"
output: 
  html_document:
    keep_md: yes
---



## Setup

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
library(ggplot2)
library(grid)
library(reshape2)
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
library(cowplot)
options(scipen = 999) #Prevent scientific notation
```



## Color palets


```r
my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442")
```





# Combine BUSCO results into tables

Combine the BUSCO proportions (percent of genes out of total in each category) for each genome together into a table using the `BUSCO_percent_to_table` bash script. Do this for each data type [genome or proteins] and BUSCO dataset [eukaryota or metazoa].


```bash
# genome vs. eukaryota_odb10
F="genome.fa.busco_eukaryota_odb10.results.txt"
./BUSCO_percent_to_table \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Montipora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Pocillopora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Porites_*/02_busco/${F} \
    ../../04_Stats_Original_Montipora_capitata_Assembly/${F} \
    genome.fa.busco_eukaryota_odb10.results.txt \
    genome_extra.fa.busco_eukaryota_odb10.results.txt \
  > combined_percent.${F}
sed -e 's@../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/@@g' \
    -e 's@/02_busco/genome.fa.busco_eukaryota_odb10.results.txt@@g' \
    -e 's@../../04_Stats_Original_Montipora_capitata_Assembly/genome.fa.busco_eukaryota_odb10.results.txt@Montipora_capitata_KBHIv1@' \
    -e 's/genome.fa.busco_eukaryota_odb10.results.txt/Montipora_capitata_KBHIv3_Chromosomes/' \
    -e 's/genome_extra.fa.busco_eukaryota_odb10.results.txt/Montipora_capitata_KBHIv3_ExtraChromosomal/' \
  -i.orig combined_percent.${F}

# genome vs. metazoa_odb10
F="genome.fa.busco_metazoa_odb10.results.txt"
./BUSCO_percent_to_table \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Montipora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Pocillopora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Porites_*/02_busco/${F} \
    ../../04_Stats_Original_Montipora_capitata_Assembly/${F} \
    genome.fa.busco_metazoa_odb10.results.txt \
    genome_extra.fa.busco_metazoa_odb10.results.txt \
  > combined_percent.${F}
sed -e 's@../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/@@g' \
    -e 's@/02_busco/genome.fa.busco_metazoa_odb10.results.txt@@g' \
    -e 's@../../04_Stats_Original_Montipora_capitata_Assembly/genome.fa.busco_metazoa_odb10.results.txt@Montipora_capitata_KBHIv1@' \
    -e 's/genome.fa.busco_metazoa_odb10.results.txt/Montipora_capitata_KBHIv3_Chromosomes/' \
    -e 's/genome_extra.fa.busco_metazoa_odb10.results.txt/Montipora_capitata_KBHIv3_ExtraChromosomal/' \
  -i.orig combined_percent.${F}

# pep vs. eukaryota_odb10
F="pep.faa.busco_eukaryota_odb10.results.txt"
./BUSCO_percent_to_table \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Montipora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Pocillopora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Porites_*/02_busco/${F} \
    ../../04_Stats_Original_Montipora_capitata_Assembly/${F} \
    pep.faa.busco_eukaryota_odb10.results.txt \
    pep_extra.faa.busco_eukaryota_odb10.results.txt \
  > combined_percent.${F}
sed -e 's@../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/@@g' \
    -e 's@/02_busco/pep.faa.busco_eukaryota_odb10.results.txt@@g' \
    -e 's@../../04_Stats_Original_Montipora_capitata_Assembly/pep.faa.busco_eukaryota_odb10.results.txt@Montipora_capitata_KBHIv1@' \
    -e 's/pep.faa.busco_eukaryota_odb10.results.txt/Montipora_capitata_KBHIv3_Chromosomes/' \
    -e 's/pep_extra.faa.busco_eukaryota_odb10.results.txt/Montipora_capitata_KBHIv3_ExtraChromosomal/' \
  -i.orig combined_percent.${F}

# pep vs. metazoa_odb10
F="pep.faa.busco_metazoa_odb10.results.txt"
./BUSCO_percent_to_table \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Montipora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Pocillopora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Porites_*/02_busco/${F} \
    ../../04_Stats_Original_Montipora_capitata_Assembly/${F} \
    pep.faa.busco_metazoa_odb10.results.txt \
    pep_extra.faa.busco_metazoa_odb10.results.txt \
  > combined_percent.${F}
sed -e 's@../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/@@g' \
    -e 's@/02_busco/pep.faa.busco_metazoa_odb10.results.txt@@g' \
    -e 's@../../04_Stats_Original_Montipora_capitata_Assembly/pep.faa.busco_metazoa_odb10.results.txt@Montipora_capitata_KBHIv1@' \
    -e 's/pep.faa.busco_metazoa_odb10.results.txt/Montipora_capitata_KBHIv3_Chromosomes/' \
    -e 's/pep_extra.faa.busco_metazoa_odb10.results.txt/Montipora_capitata_KBHIv3_ExtraChromosomal/' \
  -i.orig combined_percent.${F}
```

Combine the BUSCO counts (number of genes in each category) and percentages (of total genes in BUSCO dataset) for each genome together into a table using the `BUSCO_results_to_table` bash script. Do this for each data type [genome or proteins] and BUSCO dataset [eukaryota or metazoa].


```bash
# genome vs. eukaryota_odb10
F="genome.fa.busco_eukaryota_odb10.results.txt"
./BUSCO_results_to_table \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Montipora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Pocillopora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Porites_*/02_busco/${F} \
    ../../04_Stats_Original_Montipora_capitata_Assembly/${F} \
    genome.fa.busco_eukaryota_odb10.results.txt \
    genome_extra.fa.busco_eukaryota_odb10.results.txt \
  > combined_results.${F}
sed -e 's@../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/@@g' \
    -e 's@/02_busco/genome.fa.busco_eukaryota_odb10.results.txt@@g' \
    -e 's@../../04_Stats_Original_Montipora_capitata_Assembly/genome.fa.busco_eukaryota_odb10.results.txt@Montipora_capitata_KBHIv1@' \
    -e 's/genome.fa.busco_eukaryota_odb10.results.txt/Montipora_capitata_KBHIv3_Chromosomes/' \
    -e 's/genome_extra.fa.busco_eukaryota_odb10.results.txt/Montipora_capitata_KBHIv3_ExtraChromosomal/' \
  -i.orig combined_results.${F}

# genome vs. metazoa_odb10
F="genome.fa.busco_metazoa_odb10.results.txt"
./BUSCO_results_to_table \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Montipora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Pocillopora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Porites_*/02_busco/${F} \
    ../../04_Stats_Original_Montipora_capitata_Assembly/${F} \
    genome.fa.busco_metazoa_odb10.results.txt \
    genome_extra.fa.busco_metazoa_odb10.results.txt \
  > combined_results.${F}
sed -e 's@../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/@@g' \
    -e 's@/02_busco/genome.fa.busco_metazoa_odb10.results.txt@@g' \
    -e 's@../../04_Stats_Original_Montipora_capitata_Assembly/genome.fa.busco_metazoa_odb10.results.txt@Montipora_capitata_KBHIv1@' \
    -e 's/genome.fa.busco_metazoa_odb10.results.txt/Montipora_capitata_KBHIv3_Chromosomes/' \
    -e 's/genome_extra.fa.busco_metazoa_odb10.results.txt/Montipora_capitata_KBHIv3_ExtraChromosomal/' \
  -i.orig combined_results.${F}

# pep vs. eukaryota_odb10
F="pep.faa.busco_eukaryota_odb10.results.txt"
./BUSCO_results_to_table \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Montipora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Pocillopora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Porites_*/02_busco/${F} \
    ../../04_Stats_Original_Montipora_capitata_Assembly/${F} \
    pep.faa.busco_eukaryota_odb10.results.txt \
    pep_extra.faa.busco_eukaryota_odb10.results.txt \
  > combined_results.${F}
sed -e 's@../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/@@g' \
    -e 's@/02_busco/pep.faa.busco_eukaryota_odb10.results.txt@@g' \
    -e 's@../../04_Stats_Original_Montipora_capitata_Assembly/pep.faa.busco_eukaryota_odb10.results.txt@Montipora_capitata_KBHIv1@' \
    -e 's/pep.faa.busco_eukaryota_odb10.results.txt/Montipora_capitata_KBHIv3_Chromosomes/' \
    -e 's/pep_extra.faa.busco_eukaryota_odb10.results.txt/Montipora_capitata_KBHIv3_ExtraChromosomal/' \
  -i.orig combined_results.${F}

# pep vs. metazoa_odb10
F="pep.faa.busco_metazoa_odb10.results.txt"
./BUSCO_results_to_table \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Montipora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Pocillopora_*/02_busco/${F} \
    ../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/Porites_*/02_busco/${F} \
    ../../04_Stats_Original_Montipora_capitata_Assembly/${F} \
    pep.faa.busco_metazoa_odb10.results.txt \
    pep_extra.faa.busco_metazoa_odb10.results.txt \
  > combined_results.${F}
sed -e 's@../../../../../../0031_Coral_genomes/04_Results/2022-02-03/coral_genomes/@@g' \
    -e 's@/02_busco/pep.faa.busco_metazoa_odb10.results.txt@@g' \
    -e 's@../../04_Stats_Original_Montipora_capitata_Assembly/pep.faa.busco_metazoa_odb10.results.txt@Montipora_capitata_KBHIv1@' \
    -e 's/pep.faa.busco_metazoa_odb10.results.txt/Montipora_capitata_KBHIv3_Chromosomes/' \
    -e 's/pep_extra.faa.busco_metazoa_odb10.results.txt/Montipora_capitata_KBHIv3_ExtraChromosomal/' \
  -i.orig combined_results.${F}
```





# Load datasets into R

Load BUSCO count tables into R and melt them into a format read for ggplot.


```r
genome.metazoa <- read.table("combined_percent.genome.fa.busco_metazoa_odb10.results.txt", sep='\t', header=TRUE) %>% 
  melt(id.vars="BUSCO_categories")
genome.eukaryota <- read.table("combined_percent.genome.fa.busco_eukaryota_odb10.results.txt", sep='\t', header=TRUE) %>% 
  melt(id.vars="BUSCO_categories")
pep.metazoa <- read.table("combined_percent.pep.faa.busco_metazoa_odb10.results.txt", sep='\t', header=TRUE) %>% 
  melt(id.vars="BUSCO_categories")
pep.eukaryota <- read.table("combined_percent.pep.faa.busco_eukaryota_odb10.results.txt", sep='\t', header=TRUE) %>% 
  melt(id.vars="BUSCO_categories")
```





# Plot data

## Function to plot BUSCO results

Function to plot the percent of genes identified in each BUSCO category. Code is based on `generate_plot.py` from the BUSCO package. 


```r
plot_BUSCO_percentages = function(data, plot_title, plot_size_ratio=1, plot_family="sans")
{
  figure <- data %>%
    mutate(BUSCO_categories = factor(BUSCO_categories, c("Complete", "Single-copy", "Duplicated", "Fragmented", "Missing"))) %>%
    filter(BUSCO_categories!="Complete") %>%
    filter(variable!="Montipora_capitata_KBHIv3_Chromosomes") %>%
    filter(variable!="Montipora_capitata_KBHIv3_ExtraChromosomal") %>%
    filter(variable!="Montipora_sp1_aff_capitata_ULFMv1") %>%
    mutate(variable = factor(variable, rev(c("Montipora_capitata_KBHIv3", 
                                             "Montipora_capitata_KBHIv1", 
                                             "Montipora_capitata_WTHIv1.1", 
                                             "Montipora_cactus_SLJPv2", 
                                             "Montipora_efflorescens_SLJPv2",
                                             "Pocillopora_meandrina_KBHIv1",
                                             "Pocillopora_acuta_KBHIv2", 
                                             "Pocillopora_acuta_LBIDv1", 
                                             "Pocillopora_damicornis_SIPAv1", 
                                             "Pocillopora_verrucosa_RSSAv1",
                                             "Porites_compressa_KBHIv1",
                                             "Porites_astreoides_BBBMv1",
                                             "Porites_australiensis_SLJPv1",
                                             "Porites_lutea_OIAUv1.1",
                                             "Porites_rus_NAIDv1")))) %>%
    ggplot(aes(y = value, x = variable, fill = BUSCO_categories)) +
      geom_bar(position = position_stack(reverse = TRUE), 
               stat="identity") +
      coord_flip() +
      geom_text(aes(label = sprintf("%1.1f%%", value)), 
                size = rel(3)*plot_size_ratio, 
                fontface = "bold",
                position = position_stack(vjust = 0.5, reverse = TRUE)) +
      theme_gray(base_size = 8) +
      scale_y_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100)) +
      scale_fill_manual(values = my_colors,labels = c("Complete (C) and single-copy (S)",
                                                      "Complete (C) and duplicated (D)",
                                                      "Fragmented (F)",
                                                      "Missing (M)")) +
      ggtitle(plot_title) +
      xlab("") +
      ylab("%BUSCOs") +
      theme(plot.title = element_text(family = plot_family, 
                                      hjust=0.5, 
                                      colour = "black", 
                                      size = rel(2.2)*plot_size_ratio, 
                                      face = "bold")) +
      theme(legend.position="top",legend.title = element_blank()) +
      theme(legend.text = element_text(family = plot_family, size = rel(1.2)*plot_size_ratio)) +
      theme(legend.key.size = unit(1.5*plot_size_ratio,"line")) +
      theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) +
      theme(panel.grid.minor = element_blank()) +
      theme(panel.grid.major = element_blank()) +
      theme(axis.text.y = element_text(family = plot_family, face = "bold", size = rel(1.66)*plot_size_ratio, colour = "black")) +
      theme(axis.text.x = element_text(family = plot_family, face = "bold", size = rel(1.66)*plot_size_ratio, colour = "black")) +
      theme(axis.line = element_line(size = 1*plot_size_ratio, colour = "black")) +
      theme(axis.ticks.length = unit(0.85*plot_size_ratio, "cm")) +
      theme(axis.ticks.y = element_line(colour = "white", size = 0)) +
      theme(axis.ticks.x = element_line(colour = "#222222")) +
      theme(axis.ticks.length = unit(0.4*plot_size_ratio, "cm")) +
      theme(axis.title.x = element_text(family = plot_family, size = rel(2)*plot_size_ratio)) +
      guides(fill = guide_legend(override.aes = list(colour = NULL))) +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  return(figure)
}
```



## Plot BUSCO results for each comparison

Plot the four sets of BUSCO results that we have as separate panels (without legends), using the same color scheme as BUSCO does in its code. Extract legend from first panel and place it at the bottom of the image.


```r
plot_size_ratio <- 0.5
p1 <- plot_BUSCO_percentages(genome.metazoa, "BUSCO Metazoa dataset (Genome mode)", plot_size_ratio=plot_size_ratio)
p2 <- plot_BUSCO_percentages(genome.eukaryota, "BUSCO Eukaryota dataset (Genome mode)", plot_size_ratio=plot_size_ratio)
p3 <- plot_BUSCO_percentages(pep.metazoa, "BUSCO Metazoa dataset (Protein mode)", plot_size_ratio=plot_size_ratio)
p4 <- plot_BUSCO_percentages(pep.eukaryota, "BUSCO Eukaryota dataset (Protein mode)", plot_size_ratio=plot_size_ratio)
legend <- get_legend(p1)
grid.arrange(p1 + theme(legend.position = "none"), 
             p2 + theme(legend.position = "none"), 
             p3 + theme(legend.position = "none"), 
             p4 + theme(legend.position = "none"), 
             legend, 
            widths = c(1, 1),
            heights = c(1, 1, 0.1),
            layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)))
```

![](BUSCO_results_files/figure-html/plot_all_BUSCO_runs-1.png)<!-- -->





# Session Info


```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-apple-darwin18.7.0 (64-bit)
## Running under: macOS Mojave 10.14.6
## 
## Matrix products: default
## BLAS:   /usr/local/Cellar/openblas/0.3.18/lib/libopenblasp-r0.3.18.dylib
## LAPACK: /usr/local/Cellar/r/4.1.2/lib/R/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] cowplot_1.1.1  gridExtra_2.3  reshape2_1.4.4 ggplot2_3.3.5  dplyr_1.0.8   
## [6] readxl_1.4.0  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.8.3     highr_0.9        plyr_1.8.7       cellranger_1.1.0
##  [5] pillar_1.7.0     bslib_0.3.1      compiler_4.1.2   jquerylib_0.1.4 
##  [9] tools_4.1.2      digest_0.6.29    jsonlite_1.8.0   evaluate_0.15   
## [13] lifecycle_1.0.1  tibble_3.1.6     gtable_0.3.0     pkgconfig_2.0.3 
## [17] rlang_1.0.2      cli_3.3.0        rstudioapi_0.13  yaml_2.3.5      
## [21] xfun_0.30        fastmap_1.1.0    withr_2.5.0      stringr_1.4.0   
## [25] knitr_1.39       generics_0.1.2   vctrs_0.4.1      sass_0.4.1      
## [29] tidyselect_1.1.2 glue_1.6.2       R6_2.5.1         fansi_1.0.3     
## [33] rmarkdown_2.14   farver_2.1.0     purrr_0.3.4      magrittr_2.0.3  
## [37] scales_1.2.0     htmltools_0.5.2  ellipsis_0.3.2   colorspace_2.0-3
## [41] utf8_1.2.2       stringi_1.7.6    munsell_0.5.0    crayon_1.5.1
```

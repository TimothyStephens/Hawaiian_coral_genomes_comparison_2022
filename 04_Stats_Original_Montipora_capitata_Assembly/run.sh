

awk '$1!~"#"' Mcap.genes.gff3 | sed -e 's/;geneID=.*//' > Mcap.genes.cleaned.gff3


for F in *_odb10;
do
  ./parse_BUSCO_short_summary "${F}"/short_summary.* > "${F}".results.txt
done


scp timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/2022-01-27/02_Mcapitata_Pacuta_genomes_analysis/04_Stats_Original_Montipora_capitata_Assembly/*.results.txt .

scp timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/2022-01-27/02_Mcapitata_Pacuta_genomes_analysis/04_Stats_Original_Montipora_capitata_Assembly/*.GeneStats.tsv .





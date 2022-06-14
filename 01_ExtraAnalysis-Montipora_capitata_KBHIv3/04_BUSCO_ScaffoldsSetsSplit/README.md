# Run `BUSCO`on the *M. capitata* chromosomal and extra-chromosomal scaffolds separately. 

## Setup analysis directory

Link `fasta` files with the chromosomal and extra-chromosomal scaffolds.

```bash
ln -s ../00_databases/Montipora_capitata_KBHIv3.chromosomes.assembly.fasta
ln -s ../00_databases/Montipora_capitata_KBHIv3.extrachromosomal.assembly.fasta
```

Link to `pep` and `gff` files for predicted genes. 

```bash
ln -s ../00_databases/Montipora_capitata_KBHIv3.genes.gff3
ln -s ../00_databases/Montipora_capitata_KBHIv3.genes.cds.fna
ln -s ../00_databases/Montipora_capitata_KBHIv3.genes.pep.faa
```

Link to script for filtering text and `fasta` files.

```bash
ln -s ../../../../../02_Scripts/grepf_fasta.py
ln -s ../../../../../02_Scripts/grepf_column.py
```

Setup working environment.

```bash
conda activate py27
```

## Split predicted genes into chromosomal sets

Get the names of genes predicted on each set of chromosomes. 

Get the names of the scaffolds in each set.

```bash
grep '>' Montipora_capitata_KBHIv3.chromosomes.assembly.fasta \
  | sed -e 's/>//' > Montipora_capitata_KBHIv3.chromosomes.assembly.fasta.names
grep '>' Montipora_capitata_KBHIv3.extrachromosomal.assembly.fasta \
  | sed -e 's/>//' > Montipora_capitata_KBHIv3.extrachromosomal.assembly.fasta.names
```

Split **gff** features.

```bash
./grepf_column.py \
  -i Montipora_capitata_KBHIv3.genes.gff3 \
  -o Montipora_capitata_KBHIv3.chromosomes.genes.gff3 \
  -f Montipora_capitata_KBHIv3.chromosomes.assembly.fasta.names

./grepf_column.py \
  -i Montipora_capitata_KBHIv3.genes.gff3 \
  -o Montipora_capitata_KBHIv3.extrachromosomal.genes.gff3 \
  -f Montipora_capitata_KBHIv3.extrachromosomal.assembly.fasta.names
```

Extract gene names from each set.

```bash
cat Montipora_capitata_KBHIv3.chromosomes.genes.gff3 \
  | awk -F'\t' '$3=="transcript" {print $9}' \
  | sed -e 's/ID=//' \
  > Montipora_capitata_KBHIv3.chromosomes.genes.names

cat Montipora_capitata_KBHIv3.extrachromosomal.genes.gff3 \
  | awk -F'\t' '$3=="transcript" {print $9}' \
  | sed -e 's/ID=//' \
  > Montipora_capitata_KBHIv3.extrachromosomal.genes.names
```

Split protein `pep` file.

```bash
./grepf_fasta.py \
  -i Montipora_capitata_KBHIv3.genes.pep.faa \
  -o Montipora_capitata_KBHIv3.chromosomes.genes.pep.faa \
  -f Montipora_capitata_KBHIv3.chromosomes.genes.names

./grepf_fasta.py \
  -i Montipora_capitata_KBHIv3.genes.pep.faa \
  -o Montipora_capitata_KBHIv3.extrachromosomal.genes.pep.faa \
  -f Montipora_capitata_KBHIv3.extrachromosomal.genes.names
```

Split coding `cds` file.

```bash
./grepf_fasta.py \
  -i Montipora_capitata_KBHIv3.genes.cds.fna \
  -o Montipora_capitata_KBHIv3.chromosomes.genes.cds.fna \
  -f Montipora_capitata_KBHIv3.chromosomes.genes.names

./grepf_fasta.py \
  -i Montipora_capitata_KBHIv3.genes.cds.fna \
  -o Montipora_capitata_KBHIv3.extrachromosomal.genes.cds.fna \
  -f Montipora_capitata_KBHIv3.extrachromosomal.genes.names
```

## Run `BUSCO` analysis

Execute the `run_busco.sh` script to run `bunco` on the two scaffolds sets + associated predicted genes. Run using the `metazoa_odb10` (2021-02-24) and `eukaryota_odb10` (2020-09-10) lineages. See https://busco-data.ezlab.org/v5/data/lineages/ for the list of available lineages (version 5).

Generate BUSCO summary results file for each genome.

```bash
for F in *_odb10;
do
  ./../../../../../02_Scripts/parse_BUSCO_short_summary "${F}"/short_summary.* > "${F}".results.txt
done
```

## Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/2022-01-27/01_Hawaiian_coral_genomes/01_ExtraAnalysis-Montipora_capitata_KBHIv3/04_BUSCO_ScaffoldsSetsSplit/"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.sh*" --include="*.py" \
 --include="*.results.txt" \
 --exclude="*" \
 ${WD}/./ \
 . --dry-run
```


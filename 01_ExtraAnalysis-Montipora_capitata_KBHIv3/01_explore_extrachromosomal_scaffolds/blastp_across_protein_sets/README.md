# BLASTP across protein sets

Compare the extrachromosomal and chromosomal predicted proteins using BLASTP. Will let us see if the gene content of these scaffolds are duplicates of the main chromosomes (i.e., likely false haplotigs) or are genuine unique genes.

## Setup analysis directory

Link to each scaffold set's  `fai` file (used to extract the scaffold names in each set).

```bash
ln -s ../../00_databases/Montipora_capitata_KBHIv3.chromosomes.assembly.fasta.fai
ln -s ../../00_databases/Montipora_capitata_KBHIv3.extrachromosomal.assembly.fasta.fai
```

Link to predicted gene `gff3` file. 

```bash
ln -s ../../00_databases/Montipora_capitata_KBHIv3.genes.gff3 
```

Link to predicted gene protein files.

```bash
ln -s ../../00_databases/Montipora_capitata_KBHIv3.genes.pep.faa
ln -s ../../00_databases/Montipora_capitata_KBHIv3.genes.pep.faa.fai
```

Link to custom text parsing script for filtering table files and BLAST hits.

```bash
ln -s ../../../../../../02_Scripts/grepf_column.py
ln -s ../../../../../../02_Scripts/add_value_to_table.py
ln -s ../../../../../../02_Scripts/blast_top_hits.py
ln -s ../../../../../../02_Scripts/blast_hit_coverage.py
```

Setup bash environment.

```bash
export PATH="$PATH:/home/timothy/programs/ncbi-blast-2.10.1+/bin"
export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"

conda activate py27
```

## Extract proteins in each set

```bash
PREFIX="Montipora_capitata_KBHIv3"
OUT="extrachromosomal_vs_chromosomal_predicted_genes"
NCPUS=48
```

Extract the proteins that are predicted on either the extrachromosomal or chromosomal scaffolds into separate files.

```bash
cat "${PREFIX}.chromosomes.assembly.fasta.fai" \
  | cut -f1 \
  | ./grepf_column.py -i "${PREFIX}.genes.gff3" -f /dev/stdin \
  | awk -F'\t' '$3=="gene" {split($9,a,"="); print a[2]}' \
  | seqkit faidx --line-width 0 --infile-list /dev/stdin ${PREFIX}.genes.pep.faa \
  > "${PREFIX}.chromosomes.genes.pep.faa"
```

> 47839 proteins

```bash
cat "${PREFIX}.extrachromosomal.assembly.fasta.fai" \
  | cut -f1 \
  | ./grepf_column.py -i "${PREFIX}.genes.gff3" -f /dev/stdin \
  | awk -F'\t' '$3=="gene" {split($9,a,"="); print a[2]}' \
  | seqkit faidx --line-width 0 --infile-list /dev/stdin ${PREFIX}.genes.pep.faa \
  > "${PREFIX}.extrachromosomal.genes.pep.faa"
```

> 6545 proteins

## Align extrachromosomal proteins against chromosomal proteins

Align the extrachromosomal proteins against the chromosomal proteins using BLASTP.

```bash
makeblastdb -dbtype prot -in "${PREFIX}.chromosomes.genes.pep.faa"
blastp \
  -query "${PREFIX}.extrachromosomal.genes.pep.faa" \
  -db "${PREFIX}.chromosomes.genes.pep.faa" \
  -out "${OUT}.blastp.outfmt6" \
  -num_threads ${NCPUS} \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
```

---

Filter hits: ***e*-value < 1e-5**, **percent_ID > 75%**, and **query coverage > 75%**; take single best top hit per query sequence.

```bash
cat "${OUT}.blastp.outfmt6" \
  | awk -F'\t' '$11<1e-5 && $3>75' \
  | ./blast_hit_coverage.py -q 75 \
  | ./blast_top_hits.py -s 1 -n 1 \
  > "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit"

cut -f 1 "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit" | sort | uniq | wc -l
```

>3896 (59.53% of 6545 total)

Filter hits: ***e*-value < 1e-5**, **percent_ID > 95%**, and **query coverage > 75%**; take single best top hit per query sequence.

```bash
cat "${OUT}.blastp.outfmt6" \
  | awk -F'\t' '$11<1e-5 && $3>95' \
  | ./blast_hit_coverage.py -q 75 \
  | ./blast_top_hits.py -s 1 -n 1 \
  > "${OUT}.blastp.outfmt6.e-5.95ID.75qcov.tophit"

cut -f 1 "${OUT}.blastp.outfmt6.e-5.95ID.75qcov.tophit" | sort | uniq | wc -l
```

>2249 (34.36% of 6545 total)

Filter hits: ***e*-value < 1e-5**, **percent_ID > 95%**, and **query coverage > 95%**; take single best top hit per query sequence.

```bash
cat "${OUT}.blastp.outfmt6" \
  | awk -F'\t' '$11<1e-5 && $3>95' \
  | ./blast_hit_coverage.py -q 95 \
  | ./blast_top_hits.py -s 1 -n 1 \
  > "${OUT}.blastp.outfmt6.e-5.95ID.95qcov.tophit"

cut -f 1 "${OUT}.blastp.outfmt6.e-5.95ID.95qcov.tophit" | sort | uniq | wc -l
```

>1623 (24.80% of 6545 total)

## Parse BLASTP hits

Sort and filter the BLASTP hits so we can better see if the gene on the extrachromosomal scaffolds are subsets of the chromosomal scaffolds.

```bash
awk -F'\t' '$3=="transcript" {split($9,a,"="); print a[2]"\t"$1":"$4"-"$5":"$7}' "${PREFIX}.genes.gff3" > "${PREFIX}.genes.gff3.gene_info"
```

Reformat BLASTP results (75%/75% filtering). Output:

> qseqid   qstart-qend   qlen   sseqid   sstart-send   slen   pident   qscaffold:qstart-qstop:qstrand   sscaffold:sstart-sstop:sstrand
>
> ​                                                                                                              [            Query gene info           ]   [         Subject gene info          ]

```bash
cat "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit" \
  | awk -F'\t' '{print $1"\t"$7"-"$8"/"$13"\t"$2"\t"$9"-"$10"/"$14"\t"$3}' \
  | ./add_value_to_table.py -c 1 -a "${PREFIX}.genes.gff3.gene_info" \
  | ./add_value_to_table.py -c 3 -a "${PREFIX}.genes.gff3.gene_info" \
  > "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit.reformatted_annotated"
```

```bash
sort -k6,6 "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit.reformatted_annotated" | less
```

Reformat BLASTP results (95%/95% filtering).

```bash
cat "${OUT}.blastp.outfmt6.e-5.95ID.95qcov.tophit" \
  | awk -F'\t' '{print $1"\t"$7"-"$8"/"$13"\t"$2"\t"$9"-"$10"/"$14"\t"$3}' \
  | ./add_value_to_table.py -c 1 -a "${PREFIX}.genes.gff3.gene_info" \
  | ./add_value_to_table.py -c 3 -a "${PREFIX}.genes.gff3.gene_info" \
  > "${OUT}.blastp.outfmt6.e-5.95ID.95qcov.tophit.reformatted_annotated"
```

```bash
sort -k6,6 "${OUT}.blastp.outfmt6.e-5.95ID.95qcov.tophit.reformatted_annotated" | less
```

---

List genes along extra-chromosomal scaffolds and annotate with top hits if >95% qcov and >95% ID.

Will let us see how many genes along each scaffold have these high-similarity hits and if they are contiguous.

```bash
awk -F'\t' '$3=="transcript" {split($9,a,"="); print a[2]"\t"$1"\t"$4"\t"$5}' "${PREFIX}.genes.gff3" \
  | sort -k2,2 -k3,3n \
  | ./grepf_column.py -c 2 -f <(cut -f1 Montipora_capitata_KBHIv3.extrachromosomal.assembly.fasta.fai) \
  | ./add_value_to_table.py -a extrachromosomal_vs_chromosomal_predicted_genes.blastp.outfmt6.e-5.95ID.95qcov.tophit.reformatted_annotated \
  | less
```

---

Get the number of extra-chromosomal scaffolds with > 1 gene with a high-similarity top hit that have all of their genes (with top hits) from the same chromosomal scaffold. I.e., based on the top hit results, is the extra-chromosomal scaffold a subset of a chromosomal scaffold, or is it a chimer of multiple scaffolds. 

```bash
cut -f6,7 "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit.reformatted_annotated" \
  | sed -e 's/:[^\t]*/\t0\t1/' -e 's/:[^\t]*//' \
  | bedtools sort \
  | bedtools merge -c 4 -o count,distinct \
  > "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit.reformatted_annotated.mergedBy_query_scaffold"
```

Number of total proteins with top hits >75% qcov and >75% ID.

```bash
awk '{SUM=SUM+$4} END {print SUM}' "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit.reformatted_annotated.mergedBy_query_scaffold"
```

> 3896 (same as above; Good!)

Number of total proteins with top hits >75% qcov and >75% ID **that are on scaffolds with >1 protein with a filtered top hit**.

```bash
awk '$4>1 {SUM=SUM+$4} END {print SUM}' "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit.reformatted_annotated.mergedBy_query_scaffold"
```

> 3455 (88.68% of 3896)

Number of total proteins with top hits >75% qcov and >75% ID **that are on scaffolds with >1 protein with a filtered top hit AND that are on scaffolds with top hits to difference chromosomes** (i.e., number of genes on scaffolds that are chimers of the chromosomes).

```bash
awk '$4>1 && $5~"," {SUM=SUM+$4} END {print SUM}' "${OUT}.blastp.outfmt6.e-5.75ID.75qcov.tophit.reformatted_annotated.mergedBy_query_scaffold"
```

> 2748 (70.53% of 3896)


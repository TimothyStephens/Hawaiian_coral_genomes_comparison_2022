# BLASTN across scaffold sets

Compare the extrachromosomal and chromosomal scaffolds using BLASTN. Will let us see if these scaffolds are erroneous haplotigs or genuine sequences.

## Setup analysis directory

Link to each scaffold set  `fasta` file (used for BLASTN).

```bash
ln_loop ../../00_databases/blast/Montipora_capitata_KBHIv3.chromosomes.assembly.fasta*
ln_loop ../../00_databases/blast/Montipora_capitata_KBHIv3.extrachromosomal.assembly.fasta*
```

Link to extrachromosomal `fai` file. Used to get scaffold sizes.

```bash
ln -s ../../00_databases/Montipora_capitata_KBHIv3.extrachromosomal.assembly.fasta.fai
```

Setup bash environment.

```bash
export PATH="$PATH:/home/timothy/programs/ncbi-blast-2.10.1+/bin"
export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"

conda activate py27
```

## Align extrachromosomal against chromosomal scaffolds

```bash
PREFIX="Montipora_capitata_KBHIv3"
OUT="extrachromosomal_vs_chromosomal_scaffolds"
NCPUS=48
```

Align the extrachromosomal scaffolds against the chromosomal scaffolds.

```bash
blastn \
  -query "${PREFIX}.extrachromosomal.assembly.fasta" \
  -db "${PREFIX}.chromosomes.assembly.fasta" \
  -out "${OUT}.blastn.outfmt6" \
  -num_threads ${NCPUS} \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
```

Filter hits (*e*-value and bitscore), convert to query-based bed formatted hit regions, sort features, merge features, and get total scaffold coverage by merged hits. Will return coverage of extrachromosomal scaffolds with hits to chromosomal scaffolds.

```bash
cat "${OUT}.blastn.outfmt6" \
  | awk -F'\t' '$11<1e-5 && $12>1000 {print $1"\t"$7-1"\t"$8}' \
  | bedtools sort | bedtools merge \
  | bedtools coverage \
      -a <(awk -F'\t' '{print $1"\t0\t"$2}' "${PREFIX}.extrachromosomal.assembly.fasta.fai") \
      -b - \
  > "${OUT}.blastn.outfmt6.filtered.sortmerged.coverage.bed"
```

Output format:

>After each entry in A, reports: 
>
>1) The number of features in B that overlapped the A interval.
>2) The number of bases in A that had non-zero coverage.
>3) The length of the entry in A.
>4) The fraction of bases in A that had non-zero coverage.



```bash
zcat "${OUT}.blastn.outfmt6" \
  | awk -F'\t' '$11<1e-5 && $12>1000 && $3>95 {print $1"\t"$7-1"\t"$8}' \
  | bedtools sort | bedtools merge \
  | bedtools coverage \
      -a <(awk -F'\t' '{print $1"\t0\t"$2}' "${PREFIX}.extrachromosomal.assembly.fasta.fai") \
      -b - \
  > "${OUT}.blastn.outfmt6.filtered_95PID.sortmerged.coverage.bed"
```



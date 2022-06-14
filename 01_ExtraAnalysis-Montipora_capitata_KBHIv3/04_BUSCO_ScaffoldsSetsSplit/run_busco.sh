#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate BUSCO_v5.0.0; set -eu

NCPUS=24
PREFIX=Montipora_capitata_KBHIv3
TMPDIR="/scratch/timothy/tmp"
DB="/scratch/timothy/databases/busco"



#### Start Script
WD="$PWD"
BUSCO_TMPDIR="${TMPDIR}/busco.$(date +%s)"

run_cmd "mkdir ${BUSCO_TMPDIR}"
run_cmd "cd ${BUSCO_TMPDIR}"



## BUSCO - genome
F="genome.fa"

run_cmd "ln -s ${WD}/${PREFIX}.chromosomes.assembly.fasta ${F}"

LINEAGE="metazoa_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode genome --in ${F} --out ${F}.busco_${LINEAGE}"

LINEAGE="eukaryota_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode genome --in ${F} --out ${F}.busco_${LINEAGE}"



## BUSCO - protein
F="pep.faa"

run_cmd "ln -s $WD/${PREFIX}.chromosomes.genes.pep.faa $F"

LINEAGE="metazoa_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode protein --in ${F} --out ${F}.busco_${LINEAGE}"

LINEAGE="eukaryota_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode protein --in ${F} --out ${F}.busco_${LINEAGE}"



## BUSCO - genome
F="genome_extra.fa"

run_cmd "ln -s ${WD}/${PREFIX}.extrachromosomal.assembly.fasta ${F}"

LINEAGE="metazoa_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode genome --in ${F} --out ${F}.busco_${LINEAGE}"

LINEAGE="eukaryota_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode genome --in ${F} --out ${F}.busco_${LINEAGE}"



## BUSCO - protein
F="pep_extra.faa"

run_cmd "ln -s $WD/${PREFIX}.extrachromosomal.genes.pep.faa $F"

LINEAGE="metazoa_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode protein --in ${F} --out ${F}.busco_${LINEAGE}"

LINEAGE="eukaryota_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode protein --in ${F} --out ${F}.busco_${LINEAGE}"



## Move everything back
run_cmd "mv * ${WD}/"
run_cmd "cd ${WD}"
run_cmd "rm -fr ${BUSCO_TMPDIR}"



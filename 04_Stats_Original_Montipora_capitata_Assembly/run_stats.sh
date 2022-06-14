#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/bbmap"
export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"

#### Start Script
run_cmd "./genome_stats Mcap.genome_assembly.fa Mcap.genes.cleaned.gff3 > Mcap.GeneStats.tsv"



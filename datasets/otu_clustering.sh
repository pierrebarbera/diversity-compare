#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)

[ "$#" -eq 1 ] || die "specify dataset!"

DS=$1
DATA=${BASE}/${DS}/data

[ ! -d "${DATA}" ] && die "No such dataset."

SAMPLES=${DATA}/queries_clean
FILT_SAMPLES=${SAMPLES}/filtered

THREADS=40

OUT=${BASE}/${DS}/otu
mkdir -p ${OUT}

set -e

LOG=complete.log

cd ${OUT}
for fas in ${SAMPLES}/*.fasta; do
  name=${fas##*/}
  name=${name%.*}

  # Discard sequences containing Ns, add expected error rates
  vsearch \
      --quiet \
      --fastx_filter "${fas}" \
      --fastq_maxns 0 \
      --fastaout tmp.fasta

  # Dereplicate at the study level
  vsearch \
      --quiet \
      --derep_fulllength tmp.fasta \
      --sizeout \
      --fasta_width 0 \
      --output ${FILT_SAMPLES}/${name}.fasta

  rm tmp.fasta
done

# pool the per sample fastas
cat ${FILT_SAMPLES}/*.fasta > pooled.fasta

# global dereplication
vsearch --derep_fulllength pooled.fasta \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --output pooled.derep.fasta
rm pooled.fasta

FINAL_FASTA=pooled.derep.fasta
REPRESENTATIVES=reps.fasta

# swarming
swarm --differences 1 \
      --fastidious \
      --usearch-abundance \
      --threads ${THREADS} \
      -i ${FINAL_FASTA/.fasta/_1f.struct} \
      -s ${FINAL_FASTA/.fasta/_1f.stats} \
      -w ${REPRESENTATIVES} \
      -o ${FINAL_FASTA/.fasta/_1f.swarms} ${FINAL_FASTA}

# Sort representatives
vsearch --fasta_width 0 \
        --sortbysize ${REPRESENTATIVES} \
        --output ${FINAL_FASTA/.fasta/_1f_representatives.fas}

# rm ${REPRESENTATIVES}

#build table
${BASE}/otu_table.py --stats ${FINAL_FASTA/.fasta/_1f.stats} --swarms ${FINAL_FASTA/.fasta/_1f.swarms} ${FILT_SAMPLES}/*.fasta > otu.table

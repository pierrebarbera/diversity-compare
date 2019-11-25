#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)

set -e
shopt -s nullglob

for arx in ${BASE}/raw/*.tar; do

  name=${arx##*/}
  name=${name%.*}

  SEQDIR=${BASE}/raw/${name}
  echo "Processing ${name}..."
  cd ${SEQDIR}
  fwd=HMP.2.${name}.1.fq
  rev=HMP.2.${name}.2.fq

  # merge
  pear -f ${fwd} -r ${rev} -o ${name} --keep-original --threads 40 > pear.log

  assembl=${name}.assembled.fastq
  fasta=${name}.fasta

  # filter
  vsearch --fastq_filter ${assembl} --fastq_maxns 0 --fastaout ${fasta} --relabel_sha1 --quiet --threads 40 > vsearch.log

  # compress
  gzip -f ${fasta}

  # cleanup
  rm ${assembl}

done


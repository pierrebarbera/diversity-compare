#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)

set -e
shopt -s nullglob

cd ${BASE}

for arx in ${BASE}/raw/*.tar; do

  name=${arx##*/}
  name=${name%.*}

  FLIST="${FLIST} ${BASE}/raw/${name}/${name}.fasta.gz"

done

gappa prepare chunkify --fasta-path ${FLIST} --abundances-out-dir ${BASE}/abundance --chunks-out-dir ${BASE}/chunks --threads 40 --log-file chunkify.log --chunk-size 1000000

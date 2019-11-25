#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)

set -e
shopt -s nullglob

for arx in ${BASE}/raw/*.tar; do

  name=${arx##*/}
  name=${name%.*}

  OUT=${BASE}/raw/${name}
  mkdir -p ${OUT}

  tar -xf ${arx} --directory ${OUT}

  cd ${OUT}
  bzip2 -d *.bz2
  cd -

done


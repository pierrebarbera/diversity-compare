#!/bin/bash
BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "specify dataset!"

DS=$1
DATA=${BASE}/datasets/${DS}

[ ! -d "${DATA}" ] && die "No such dataset."

WORKDIR=${BASE}/workdir/${DS}/scrapp_boot_1000

JPDIR=${DATA}/place/samples
SEQDIR=${JPDIR}

mkdir -p ${WORKDIR}

# set -e
shopt -s nullglob

N=10
(
for jplace in ${JPDIR}/*.jplace; do

  name=${jplace##*/}
  name=${name%.*}

  OUT=${WORKDIR}/${name}
  QRY=${SEQDIR}/${name}.fasta.gz

  ((i=i%N)); ((i++==0)) && wait
  ~/scrapp/scrapp.py --cluster-above 1000 --bootstrap --jplace ${jplace} --alignment ${QRY} --work-dir ${OUT} --parallel threads --num-threads 4 --min-weight 0.5 &

done
)

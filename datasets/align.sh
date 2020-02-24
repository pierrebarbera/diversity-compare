#!/bin/bash

THREADS=4

die () {
    echo >&2 "$@"
    exit 1
}

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)

[ "$#" -eq 1 ] || die "specify dataset!"

DS=$1
DSBASE=${BASE}/${DS}
DATA=${DSBASE}/data
SAMPLES=${DATA}/samples

[ ! -d "${DATA}" ] && die "No such dataset."

TREEDIR=${DSBASE}/tree
TREE=${TREEDIR}/reference.raxml.bestTree
ALI=${DATA}/reference.phylip

[ ! -f "${TREE}" ] && die "Reference tree not found: ${TREE}"

QUERIES=${DATA}/chunks/chunk_0.fasta

[ ! -f "${QUERIES}" ] && die "Queries not found: ${QUERIES}"

ALIGND=${DSBASE}/align

mkdir -p ${ALIGND}
cd ${ALIGND}

papara -t ${TREE} -s ${ALI} -q ${QUERIES} -j ${NUM_TASKS} -r -n chunk_0

# split results for epa
epa-ng --split ${ALI} papara_alignment.chunk_0
gzip query.fasta
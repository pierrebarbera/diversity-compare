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

[ ! -d "${DATA}" ] && die "No such dataset."

TREEDIR=${DSBASE}/tree
TREE=${TREEDIR}/reference.raxml.bestTree
ALIGND=${DSBASE}/align

[ ! -f "${TREE}" ] && die "Reference tree not found: ${TREE}"

QRY=${ALIGND}/query.fasta

[ ! -f "${QRY}" ] && die "Queries not found: ${QRY}"

REF=${ALIGND}/reference.fasta
MODEL=${TREEDIR}/reference.raxml.bestModel

PLACEDIR=${DSBASE}/place

mkdir -p ${PLACEDIR}

epa-ng --ref-msa ${REF} --tree ${TREE} --query ${QRY} --model ${MODEL} --out-dir ${PLACEDIR}

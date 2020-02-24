#!/bin/bash

THREADS=4

die () {
    echo >&2 "$@"
    exit 1
}

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)

[ "$#" -eq 1 ] || die "specify dataset!"

DS=$1
DATA=${BASE}/${DS}/data

[ ! -d "${DATA}" ] && die "No such dataset."

mkdir -p ${DATA}/chunks
mkdir -p ${DATA}/maps

gappa prepare chunkify --fasta-path ${DATA}/samples --chunks-out-dir ${DATA}/chunks --abundances-out-dir ${DATA}/maps --threads $THREADS --chunk-size 10000000

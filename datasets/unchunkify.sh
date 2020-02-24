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


PLACEDIR=${DSBASE}/place

[ ! -d "${PLACEDIR}" ] && die "No place directory: maybe do placement first? Consult README.md!"

mkdir -p ${PLACEDIR}/samples

gappa prepare unchunkify --abundances-path ${DATA}/maps --jplace-path ${PLACEDIR}/epa_result.jplace --out-dir ${PLACEDIR}/samples

#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)

DATA=${BASE}/../data

set -e

cd ${BASE}

pbgappa prepare unchunkify --sequence-path ${BASE}/../query/query.fasta.gz --jplace-path ${BASE}/epa_result.jplace --abundances-path ${DATA}/abundance --out-dir ${BASE}/samples --threads 40 --log-file unchunkify.log

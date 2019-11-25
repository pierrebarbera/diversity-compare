#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)

cd $BASE
wget --recursive --no-parent --no-directories -P raw ftp://public:hmp2_ftp@ftp.broadinstitute.org/raw/HMP2/16S/2018-02-06

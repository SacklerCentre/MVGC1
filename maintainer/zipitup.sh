#!/bin/sh

# run from MVGC root directory!

# $1 is version (eq 1.0), $2 is directory (eg /tmp)

cd .. && zip -r -v $2/mvgc_v$1.zip mvgc_v$1 -x mvgc_v$1/testing\* mvgc_v$1/deprecated\* mvgc_v$1/maintainer\*

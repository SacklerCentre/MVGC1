#!/bin/sh

# run from MVGC root directory!

# $1 is version (eq 1.0), $2 is directory (eg /tmp)

cd .. && tar czvf $2/mvgc_v$1.tar.gz --exclude=testing --exclude=deprecated --exclude=maintainer mvgc_v$1/*

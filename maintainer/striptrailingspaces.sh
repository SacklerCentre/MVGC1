#!/bin/sh

# run e.g. as
#
# striptrailingspaces.sh \*.m

for f in $1
do
    echo "Processing $f"
    sed 's/[ \t]*$//' $f > $f.tmp && mv $f.tmp $f
done

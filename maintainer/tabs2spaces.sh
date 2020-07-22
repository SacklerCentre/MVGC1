#!/bin/sh

# run e.g. as
#
# tabs2spaces.sh \*.m

for f in $1
do
    echo "Processing $f"
    expand -t 4 $f > $f.tmp && mv $f.tmp $f
done

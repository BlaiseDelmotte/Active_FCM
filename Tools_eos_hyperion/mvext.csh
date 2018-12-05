#!/bin/sh
for ID in `ls *$1`; do
mv $ID `basename $ID $1`$2;
done

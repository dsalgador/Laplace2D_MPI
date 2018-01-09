#!/bin/bash
infile=$1
outfile=$2
awk 'BEGIN{i=0}{i++;if (i%2==0 || i==1) print $1}' $infile > $outfile
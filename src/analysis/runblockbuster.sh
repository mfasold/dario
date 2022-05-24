#!/bin/bash

timeout 10m $1analysis/blockbuster.x -scale 0.4 $2unknown.reads > $2unknown.clusters

timeout 10m $1analysis/blockbuster.x -scale 0.4 $2ncRNAs.reads > $2ncRNAs.clusters

$1analysis/flagKnownClusters.pl -c $2ncRNAs.clusters -a $3 -p 1 > $2ncRNAs.clusters.flagged
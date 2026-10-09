#!/bin/sh
f=_sources/gaia/$1.tsv; [ -s $f ] && exit 0
curl -s -m 60 "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/355/gaiadr3&-c=$2+$3&-c.rs=30&-out=Source,Gmag,RUWE&-out.add=_r&-sort=_r" > $f

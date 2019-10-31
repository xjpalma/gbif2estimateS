#!/usr/bin/bash

#source activate py3.7

python ../dist/gbif2estimateS.py -v -i occurrence.txt -o output/output.dat -geojson gis/domain.geojson

#source deactivate

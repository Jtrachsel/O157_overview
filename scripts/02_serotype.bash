#!/bin/bash


# unzip downloaded genomes
# run ectyper on them
# copy results to final location
# cleanup


mkdir -p ./data/genomes/tmp_ectyper



parallel -j 20 'gunzip {}' ::: ./data/genomes/*gz


# this generates the file "output.tsv"
ectyper -i ./data/genomes/ -o ./data/genomes/tmp_ectyper -c 20


# remove the headers from the file and append to all_sero file
sed 1d ./data/genomes/tmp_ectyper/output.tsv >> ./data/all_sero.tsv


rm -r ./data/genomes/tmp_ectyper



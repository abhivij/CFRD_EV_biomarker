#!/bin/sh

cd /srv/scratch/vafaeelab/AbhishekVijayan/exoCFRD_miRNA_2023/UMA12487
mkdir UMA12487_single_level

for f_name in UMA12487/*/*
do
    cp $f_name UMA12487_single_level/
done
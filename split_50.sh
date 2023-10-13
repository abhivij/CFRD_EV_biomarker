#!/bin/sh

cd /srv/scratch/vafaeelab/AbhishekVijayan/exoCFRD_miRNA/WAT10166-WAT10180-WAT9739-single-level/
i=0
for f_name in *
do
 dir_name="dir"$(( (i/50) + 1 ))
 echo $dir_name
 echo $f_name
 dir_path="../"$dir_name
 if [ ! -d "$dir_path" ]; then
  mkdir $dir_path
 fi
 cp $f_name $dir_path"/"
 i=$((i+1))
done

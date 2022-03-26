#!/bin/bash

p=$1

wget "https://www.encodeproject.org/files/${p}/@@download/${p}.bam" -O /oak/stanford/groups/engreitz/Users/kmualim/ABC_data/${p}.bam
samtools view -F 780 -q 30 -u /oak/stanford/groups/engreitz/Users/kmualim/ABC_data/${p}.bam | samtools sort /dev/stdin -o /oak/stanford/groups/engreitz/Users/kmualim/ABC_data/$p.se.bam -@ 20 
samtools index /oak/stanford/groups/engreitz/Users/kmualim/ABC_data/$p.se.bam

#!/bin/bash

module load singularity

cd /projects/carter-lab/hepatocytes/ChromHMM/data/bamfiles

targetdir=$/projects/carter-lab/hepatocytes/ChromHMM/data/bedfiles


for x in *.sorted.B6co.bam ; do
    echo "print current:$x";
    #echo "$targetdir${x%.bam}.bed";
    bedtools bamToBed -i "$x" > "$targetdir${x%.bam}.bed";
done
echo "done"

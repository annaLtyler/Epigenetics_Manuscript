#!/bin/bash

module load singularity

inputdir=$/projects/carter-lab/hepatocytes/ChromHMM/data/bamfiles
targetdir=$/projects/carter-lab/hepatocytes/ChromHMM/data/bedfiles

for x in *.sorted.B6co.bam ; do
    echo "print current:$x";
    #echo "$targetdir${x%.bam}.bed";
    singularity exec bedtools_v2.27.1dfsg-4-deb_cv1.sif 'bedtools bamToBed -i' "$x" > "$targetdir${x%.bam}.bed";
done
echo "done"

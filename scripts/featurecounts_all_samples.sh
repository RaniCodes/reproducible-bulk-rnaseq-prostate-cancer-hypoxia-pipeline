#!/bin/bash

# Go to your aligned reads folder
cd /path_to_folder_where_bamfiles_are_stored/alignedreads

# Loop over all BAM files
for bam in *.bam; do
    start=$(date +%s)  # start time

    echo "Processing $bam ..."
    featureCounts -s 0 -a /path_to_annotationfile/Homo_sapiens.GRCh38.114.gtf \
        -o /path_to_quants_folder/quants/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes

    echo "Completed $bam in $runtime minutes."
    echo "------------------------------------"
done

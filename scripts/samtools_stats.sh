#!/usr/bin/bash

# Assign list of SAM files
FILES=output/part_1/*.sam

for file in $FILES
do
	# Get filename (with path) but the extension
	fname="${file%.*}"

	# Convert SAM file to BAM
	samtools sort "$file" > "$fname".bam

	# Index the BAM file
	samtools index "$fname".bam

	# Get alignment statistics
	samtools stats "$fname".bam > "$fname"_stats.txt

	# Get flagstat statistics
	samtools flagstat "$fname".bam > "$fname"_flagstat.txt
done

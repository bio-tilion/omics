---
title: "Omics coursework - Part 1"
author: "Alberto Locca"
output: html_document
---

## Working environment
I have performed both the NGS practical and all the coursework locally on my machine, therefore I had prepared a new `conda` environment and installed the necessary software with the following command:
```
conda install cutadapt bowtie2 samtools fastqc multiqc igv
```

I then saved the environment for future use:
```
conda list --explicit > omics-env.txt
```

Calling the following command will recreate my environment with all the necessary software and its dependencies:
```
conda env create --file omics-env.txt
```
   

## Original reads from the practical
The FastQC report on the `Negative.fq` from the practical only flagged `Per base sequence quality` and `Per base N content`.

![Per base sequence quality](../input/part_1/fastq/trimmed_Negative_fastqc/Images/per_base_quality.png)

The graph clearly shows that the reads all have very high quality in all positions but the first five and last four positions.

![Per base N content](../input/part_1/fastq/trimmed_Negative_fastqc/Images/per_base_n_content.png)

The N content graph shows that those above mentioned positions are N calls, which are present in all reads. N calls occur when during sequencing there is not enough confidence to make a base call (A, T, C, or G), so a "wildcard" base is introduced.

## Mapping options
During the practical, we have always aligned the reads with the `--end-to-end` setting, which tries to align each read from the first to the last position to the reference genome. In this case, it failed to provide any alignment because it could not find any match to the Ns starting sequences.

One possible alternative is to perform the alignment using the `--local` setting, which performs a local alignment resulting in possible clipping of the extremities. In this case it should be able to provide an alignment and *soft-clip* the polyN head and tail.

```
bowtie2 --local --all -x input/part_1/genome/AFPN02.1_merge -q input/part_1/fastq/trimmed_Negative.fq -S output/part_1/Negative_local.sam >& output/part_1/Negative_local_bowtie_stats.txt
```

Unfortunately, the output of the bowtie2 command (`Negative_local_bowtie_stats.txt`) confirms that none of the reads have been aligned:
```
1076320 reads; of these:
  1076320 (100.00%) were unpaired; of these:
    1076320 (100.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
0.00% overall alignment rate
```

Another possible solution is to pass to `bowtie2` the number of bases to trim, using `--trim5` (or `-5`) and `--trim3` (or `-3`) followed by the number of positions. In this case, I have kept the original alignment setting and set to remove the N bases reported by the FastQC report -- the first 5 bases, and the last 4 of each read. 

```
bowtie2 --end-to-end --all -x input/part_1/genome/AFPN02.1_merge -q input/part_1/fastq/trimmed_Negative.fq -5 5 -3 4 -S output/part_1/Negative_trim.sam >& output/part_1/Negative_trim_bowtie_stats.txt
```

This time, the output (`Negative_trim_bowtie_stats.txt`) confirms that it has produced an alignment to the reference genome with an overall alignment rate of 99.97%:
```
1076320 reads; of these:
  1076320 (100.00%) were unpaired; of these:
    308 (0.03%) aligned 0 times
    1020237 (94.79%) aligned exactly 1 time
    55775 (5.18%) aligned >1 times
99.97% overall alignment rate
```

This is confirmation that the alignment failed because of the N calls in the reads.

As a slightly different alternative, I have also tried to include trimming of the N calls into the pipeline. This could be achieved during the "de-multiplexing" step using `cutadapt`, by passing the option `--trim-n`, which removes Ns at both ends of the reads.

In this case, I have used the `trimmed_Negative.fq` file which was already the output of `cutadapt` in the practical, but the trim option could have been combined with the de-multiplexing command too.
```
cutadapt --trim-n -o input/part_1/fastq/trimmed_N_Negative.fq input/part_1/fastq/trimmed_Negative.fq
```

A new FastQC report now shows that all the reads have sufficient quality in all positions. 

![Per base sequence quality after N trimming with `cutadapt`](../input/part_1/fastq/trimmed_N_Negative_fastqc/Images/per_base_quality.png)

The following is the `summary.txt` file confirming that all tests have passed:
```
PASS    Basic Statistics        trimmed_N_Negative.fq
PASS    Per base sequence quality       trimmed_N_Negative.fq
PASS    Per sequence quality scores     trimmed_N_Negative.fq
PASS    Per base sequence content       trimmed_N_Negative.fq
PASS    Per sequence GC content trimmed_N_Negative.fq
PASS    Per base N content      trimmed_N_Negative.fq
PASS    Sequence Length Distribution    trimmed_N_Negative.fq
PASS    Sequence Duplication Levels     trimmed_N_Negative.fq
PASS    Overrepresented sequences       trimmed_N_Negative.fq
PASS    Adapter Content trimmed_N_Negative.fq
```

I have then produced a new alignment with the following code:
```
bowtie2 --end-to-end --all -x input/part_1/genome/AFPN02.1_merge -q input/part_1/fastq/trimmed_N_Negative.fq -S output/part_1/Negative_cutadapt.sam >& output/part_1/Negative_cutadapt_bowtie_stats.txt
```

The output (`Negative_cutadapt_bowtie_stats.txt`) reports the same results as the alignment achieved with `bowtie2` and the trim options:
```
1076320 reads; of these:
  1076320 (100.00%) were unpaired; of these:
    308 (0.03%) aligned 0 times
    1020237 (94.79%) aligned exactly 1 time
    55775 (5.18%) aligned >1 times
99.97% overall alignment rate
```

## Final mapping statistics
In order to obtain all the mapping statistics, I wrote this short script -- `script/samtools_stats.sh` -- that iterates through the SAM alignment files:
```
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
```

Finally, I have produced a summary of all statistics with the following command:
```
multiqc input/part_1/ output/part_1/ -f -o output/part_1/
```

The report can be accessed at `output/part_1/multiqc_report.html`.

## Conclusion
![IGV snapshot of aligned reads](../output/part_1/igv_snapshot.png)
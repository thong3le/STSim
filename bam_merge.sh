#!/bin/bash

B1=subclone_at_leave_a_10x_sorted.bam
B2=subclone_at_leave_b_10x_sorted.bam
B3=subclone_at_leave_c_10x_sorted.bam
B4=subclone_at_leave_d_10x_sorted.bam
B5=subclone_at_leave_e_10x_sorted.bam

OUT=somatic_40x.bam

# create header for bam files that will be merged
# the following header is for merging B2, B3, B4, B5
samtools view -H $B2 > header.sam
samtools view -H $B3 | grep @RG >> header.sam
samtools view -H $B4 | grep @RG >> header.sam
samtools view -H $B5 | grep @RG >> header.sam

# merge the bam files

samtools merge -h header.sam $OUT $B2 $B3 $B4 $B5

samtools index $OUT 

rm header.sam

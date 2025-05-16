#!/usr/bin/env nextflow


/*
 * Given an input bed file, separate the content into those that intersect with
 * the records in another bed file and those that do not intersect.
 */


process separationByOverlap {
    label 'bedtools'
    publishDir params.publish_dir, mode: 'copy'

    cpus 1
    memory 10.GB
    input:
        path input_bed
        path separation_bed
    
    output:
        path "${input_bed}.intersected"
        path "${input_bed}.not_intersected"

    script:
    """
    bedtools intersect -a $input_bed -b $separation_bed -u -f 0.9 > ${input_bed}.intersected
    bedtools intersect -a $input_bed -b $separation_bed -v -f 0.9 > ${input_bed}.not_intersected
    """
}
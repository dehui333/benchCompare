#!/usr/bin/env nextflow

/*
 * Given a list of read ids and a bam for these reads, output the alignment start 
 * and end of of each read in the list. (taking primary alignment)
 */


process getReadIds {
    label 'samtools'

    cpus 1
    memory 10.GB
    input:
        path bef_bam, name: 'bef.bam' 
        path aft_bam, name: 'aft.bam'
    
    output:
        path 'ids.txt'

    script: 
    """
    samtools view $bef_bam | cut -f1 > duplicate_ids.txt
    samtools view $aft_bam | cut -f1 >> duplicate_ids.txt
    sort duplicate_ids.txt | uniq > ids.txt
    """
}

process filterBamRegion {
    label 'samtools'

    cpus 1
    memory 10.GB
    input:
        path bam
        path idx 
        path bed
    
    output:
        path 'filtered_region.bam', emit: bam

    script: 
    """
    samtools view -b --region-file $bed -F 3840 -o filtered_region.bam $bam  
    """
}


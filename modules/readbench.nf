#!/usr/bin/env nextflow

/*
 * Run q100bench reads evaluation.
 * 
 */

process readbench {
    label 'readbench'

    cpus 1
    memory 4.GB
    input:
        path activate
        path config
        path ref 
        tuple val(idx), path(bam), path(bam_index)

    output:
        tuple val(idx), path('out/readbench.readerrors.txt'), emit: errors
        tuple val(idx), path('out/readbench.generalstats.txt'), emit: stats
        

    script: 
    """
    source $activate
    readbench -b $bam -r $ref -p out -B q100bench -R readbench -c $config --nomononucs
    """
}

process readbenchWithBed {
    label 'readbench'

    cpus 1
    memory 4.GB
    input:
        path activate
        path config
        path ref 
        tuple val(idx), path(bam), path(bam_index)
        path bed

    output:
        tuple val(idx), path('out/readbench.readerrors.txt'), emit: errors
        tuple val(idx), path('out/readbench.generalstats.txt'), emit: stats
        

    script: 
    """
    source $activate
    readbench -b $bam -r $ref --regions $bed -p out -B q100bench -R readbench -c $config --nomononucs
    """
}


/*
 * Filter bam to contain only specific reads. 
 * 
 */

process filterBam {
    conda '/home/lindehui/miniforge3/envs/eval'
    cpus 4
    memory 4.GB
    //publishDir 'filterBam_result', mode: 'copy'
    input:
        tuple val(idx), path(ids), path(bam)
    
    output:
        tuple val(idx), path('filtered.bam'), path('filtered.bam.bai')
        

    script: 
    """
    samtools view -N $ids -F 3840 $bam -@$task.cpus -b -o filtered.bam
    samtools index -@$task.cpus filtered.bam
    """
}

workflow {
    //readbench(params.activate, params.config, params.bam, "${params.bam}.bai", params.ref)
    //filterBam(params.ids, params.bam)
}
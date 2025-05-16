#!/usr/bin/env nextflow
include { filterBam as filterBamIdsBefore } from './modules/readbench.nf'
include { filterBam as filterBamIdsAfter } from './modules/readbench.nf'
include { readbench as readbenchBefore } from './modules/readbench.nf'
include { readbench as readbenchAfter } from './modules/readbench.nf'
include { readbenchWithBed as readbenchWithBedBefore } from './modules/readbench.nf'
include { readbenchWithBed as readbenchWithBedAfter } from './modules/readbench.nf'
include { detectMovedWorkflow } from './modules/detectMoved.nf'
include { checkOutput } from './modules/checkOutput.nf'
include { splitWorkloadWorkflow } from './modules/splitWorkload.nf'
include { compareBench } from './modules/compareBench.nf'
include { mergeOutputs } from './modules/mergeOutputs.nf'
include { filterBamRegion as filterBamRegionBefore } from './modules/utilities.nf'
include { filterBamRegion as filterBamRegionAfter } from './modules/utilities.nf'
include { getReadIds } from './modules/utilities.nf'
include { separationByOverlap as separateAdded } from './modules/separationByOverlap.nf'
include { separationByOverlap as separateRemaining } from './modules/separationByOverlap.nf'
include { getStats as getBeforeStats } from './modules/mergeOutputs.nf'
include { getStats as getAfterStats } from './modules/mergeOutputs.nf'
include { splitMovedByType } from './modules/splitMovedByType.nf'
include { combineWithTriobin } from './modules/splitMovedByType.nf'

workflow {
    bef_bam = params.bef_bam
    aft_bam = params.aft_bam

    if (params.include_bed != null) {
        filterBamRegionBefore(params.bef_bam, "${params.bef_bam}.bai", params.include_bed)
        filterBamRegionAfter(params.aft_bam, "${params.aft_bam}.bai", params.include_bed)
        bef_bam = filterBamRegionBefore.out.bam
        aft_bam = filterBamRegionAfter.out.bam
    }

    getReadIds(bef_bam, aft_bam)

    splitWorkloadWorkflow(bef_bam, aft_bam, getReadIds.out, params.split_size)
    idx_ch = splitWorkloadWorkflow.out.bef_bam_ch.count().map { n -> (0..n-1) }.flatten()
    bef_bam_ch = idx_ch.merge(splitWorkloadWorkflow.out.bef_bam_ch)
    aft_bam_ch = idx_ch.merge(splitWorkloadWorkflow.out.aft_bam_ch)
    ids_ch = idx_ch.merge(splitWorkloadWorkflow.out.ids_ch)
    

    detectMovedWorkflow(bef_bam_ch, aft_bam_ch, ids_ch)
    
    
    filterBamIdsBefore(detectMovedWorkflow.out.not_moved_ids.join(bef_bam_ch))
    filterBamIdsAfter(detectMovedWorkflow.out.not_moved_ids.join(aft_bam_ch))

    
    if (params.q100bench_bed != null) {
        readbenchWithBedBefore(params.q100bench_activate, params.q100bench_config, params.ref, filterBamIdsBefore.out, params.q100bench_bed)
        readbenchWithBedAfter(params.q100bench_activate, params.q100bench_config, params.ref, filterBamIdsAfter.out, params.q100bench_bed)
        benchBeforeOut = readbenchWithBedBefore.out
        benchAfterOut = readbenchWithBedAfter.out
    } else {
        readbenchBefore(params.q100bench_activate, params.q100bench_config, params.ref, filterBamIdsBefore.out)
        readbenchAfter(params.q100bench_activate, params.q100bench_config, params.ref, filterBamIdsAfter.out)
        benchBeforeOut = readbenchBefore.out
        benchAfterOut = readbenchAfter.out
    }
    
    
    
    compareBench(benchBeforeOut.errors.join(benchAfterOut.errors))

    checkOutput(compareBench.out.join(benchBeforeOut.errors).join(benchAfterOut.errors))


    remaining_errors = compareBench.out.map { tuple -> tuple[1] }.collect()
    removed_errors = compareBench.out.map { tuple -> tuple[2] }.collect()
    added_errors = compareBench.out.map { tuple -> tuple[3] }.collect()
    moved = detectMovedWorkflow.out.moved.map { tuple -> tuple[1] }.collect()
    not_moved = detectMovedWorkflow.out.not_moved.map { tuple -> tuple[1] }.collect()
    unmapped_ids = detectMovedWorkflow.out.unmapped_ids.map { tuple -> tuple[1] }.collect()
    not_in_bam_ids = detectMovedWorkflow.out.not_in_bam_ids.map { tuple -> tuple[1] }.collect()
    debug = checkOutput.out.map { tuple -> tuple[1] }.collect()
    
    mergeOutputs(remaining_errors,
     removed_errors,
     added_errors,
     moved,
     not_moved,
     unmapped_ids,
     not_in_bam_ids,
     debug)
    if (params.separation_bed != null) {
        separateAdded(mergeOutputs.out.added_errors, params.separation_bed)
        separateRemaining(mergeOutputs.out.remaining_errors, params.separation_bed)
    }
    

    before_stats = benchBeforeOut.stats.map { tuple -> tuple[1] }.collect()
    after_stats = benchAfterOut.stats.map { tuple -> tuple[1] }.collect()
    getBeforeStats(before_stats, 'before')
    getAfterStats(after_stats, 'after')

    splitMovedByType(mergeOutputs.out.moved)

    if (params.bef_triobin != null && params.aft_triobin != null) {
        combineWithTriobin(splitMovedByType.out.hap_flip, params.bef_triobin, params.aft_triobin)
    }
    
}
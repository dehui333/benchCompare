#!/usr/bin/env nextflow
include { getAlignmentSpan as getBeforeAlignmentSpan} from './getAlignmentSpan.nf'
include { getAlignmentSpan as getAfterAlignmentSpan} from './getAlignmentSpan.nf'

/*
 * Given alignment spans before and after, find out which of the
 * reads changed their alignment.
 */

process detectMoved {

    cpus 1
    memory 4.GB

    input:
        tuple val(idx), path('before_span.tsv'), path('after_span.tsv'), path('ids.txt')
    
    output:
        tuple val(idx), path('moved.tsv'), emit: moved
        tuple val(idx), path('not_moved.tsv'), emit: not_moved
        tuple val(idx), path('not_moved.ids'), emit: not_moved_ids
        tuple val(idx), path('unmapped_in_both.ids'), emit: unmapped_ids
        tuple val(idx), path('not_in_bam.ids'), emit: not_in_bam_ids
    script: 
    """
#!/usr/bin/env python
PROP_OVERLAP = 0.95
def construct_dictionary(span_tsv:str):
    dic = {}
    with open(span_tsv, 'r') as f:
        for line in f:
            parts = line.strip().split()
            read_id = parts[3]
            ref_id = parts[0]
            if parts[1] == '*':
                ref_start = parts[1]
                ref_end = parts[2]
            else:
                ref_start = int(parts[1])
                ref_end = int(parts[2])
            dic[read_id] = (ref_id, ref_start, ref_end)
    return dic
with open('ids.txt', 'r') as f:
    all_ids = {line.strip() for line in f}

before_dict = construct_dictionary('before_span.tsv')
after_dict = construct_dictionary('after_span.tsv')

before_ids = set(before_dict.keys())
after_ids = set(after_dict.keys())

not_in_both_ids = all_ids - (before_ids | after_ids)
only_in_before_ids = before_ids - after_ids
only_in_after_ids = after_ids - before_ids

in_both_ids = before_ids & after_ids

with open('not_in_bam.ids', 'w') as f:
    for read_id in only_in_before_ids:
        f.write(f'{read_id},ONLY_IN_BEFORE\\n')
    for read_id in only_in_after_ids:
        f.write(f'{read_id},ONLY_IN_AFTER\\n')
    for read_id in not_in_both_ids:
        f.write(f'{read_id},NOT_IN_BOTH\\n')
    


with open('moved.tsv', 'w') as moved, open('not_moved.tsv', 'w') as not_moved, \
 open('not_moved.ids', 'w') as not_moved_ids, open('unmapped_in_both.ids', 'w') as unmapped_ids:
    for read_id in in_both_ids:
        before_ref_id, before_ref_start, before_ref_end = before_dict[read_id]
        after_ref_id, after_ref_start, after_ref_end = after_dict[read_id]
        if before_ref_id == 'UNMAPPED' and after_ref_id == 'UNMAPPED':
            unmapped_ids.write(f'{read_id}\\n')
            continue

        if before_ref_id != after_ref_id:
            moved.write(f'{before_ref_id}\\t{before_ref_start}\\t{before_ref_end}\\t' +
            f'{after_ref_id}\\t{after_ref_start}\\t{after_ref_end}\\t{read_id}\\n')
        else:
            after_alignment_len = after_ref_end - after_ref_start
            len_overlap_with_bef_span = min(before_ref_end, after_ref_end) - max(before_ref_start, after_ref_start)
            if len_overlap_with_bef_span >= PROP_OVERLAP * after_alignment_len :
                not_moved.write(f'{before_ref_id}\\t{before_ref_start}\\t{before_ref_end}\\t' +
                    f'{after_ref_id}\\t{after_ref_start}\\t{after_ref_end}\\t{read_id}\\n')
                not_moved_ids.write(f'{read_id}\\n')
            else:
                moved.write(f'{before_ref_id}\\t{before_ref_start}\\t{before_ref_end}\\t' +
                    f'{after_ref_id}\\t{after_ref_start}\\t{after_ref_end}\\t{read_id}\\n')
    """
}

workflow detectMovedWorkflow {
    take:
        bef_bam_ch
        aft_bam_ch
        ids_ch
    main:

        getBeforeAlignmentSpan(bef_bam_ch.join(ids_ch))
        getAfterAlignmentSpan(aft_bam_ch.join(ids_ch))

        detectMoved(getBeforeAlignmentSpan.out.join(getAfterAlignmentSpan.out).join(ids_ch))
    emit:
        moved = detectMoved.out.moved
        not_moved = detectMoved.out.not_moved
        not_moved_ids = detectMoved.out.not_moved_ids
        unmapped_ids = detectMoved.out.unmapped_ids
        not_in_bam_ids = detectMoved.out.not_in_bam_ids

}
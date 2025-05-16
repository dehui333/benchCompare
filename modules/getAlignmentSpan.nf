#!/usr/bin/env nextflow

/*
 * Given a list of read ids and a bam for these reads, output the alignment start 
 * and end of of each read in the list. (taking primary alignment)
 */


process getAlignmentSpan {
    label 'pysam'

    cpus 1
    memory 4.GB
    input:
        tuple val(idx), path(bam), path(ids)
    
    output:
        tuple val(idx), path('alignment_spans.tsv')

    script: 
    """
#!/usr/bin/env python
import pysam

with open('$ids', 'r') as ids_file:
    set_of_ids = {line.strip() for line in ids_file}


with pysam.AlignmentFile('$bam') as bam, open('alignment_spans.tsv', 'w') as out:
    for record in bam.fetch(until_eof=True):
        if record.query_name not in set_of_ids:
            continue
        if record.is_secondary or record.is_supplementary:
            continue

        query_name = record.query_name
        if record.is_unmapped:
            reference_name = 'UNMAPPED'
            reference_start = '*'
            reference_end = '*'
        else:
            reference_name = record.reference_name
            reference_start = record.reference_start
            reference_end = record.reference_end

        out.write(f'{reference_name}\\t{reference_start}\\t{reference_end}\\t{query_name}\\n')
    """
}
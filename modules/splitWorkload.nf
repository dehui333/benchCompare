#!/usr/bin/env nextflow


/*
 * Given pairs of bef/aft bams, and a list of read ids, split into smaller   
 * size (bef_bam, aft_bam, ids) files.
 */
// TODO: do not really split if size of input is <= split size
process splitWorkload {
    label 'pysam'
    cpus 1

    input:
    path bef_bam, name: 'bef.bam'
    path aft_bam, name: 'aft.bam'
    path ids
    val size

    output:
    path 'bef_bam_*', emit: bef_bam
    path 'aft_bam_*', emit: aft_bam
    path 'ids_*', emit: ids

    script:
    """
    #!/usr/bin/env python
    import pysam

    def split_ids(path:str) -> dict[str, int]:
        dic = {}
        current_output_number = 0
        current_output_file = open(f'ids_{current_output_number}.txt', 'w')

        with open('$ids', 'r') as f:
            for i, line in enumerate(f):
                read_id = line.strip()

                number = i // $size
                if number != current_output_number:
                    current_output_file.close()
                    current_output_number = number
                    current_output_file = open(f'ids_{current_output_number}.txt', 'w')
                
                current_output_file.write(line)
                dic[read_id] = number
        current_output_file.close()
        return dic
    
    def split_bam(path:str, out_prefix:str, rid_to_number:dict[str, int]) -> None:
        
        dic_of_files = {}
        with pysam.AlignmentFile(path) as bam:
            for record in bam.fetch(until_eof=True):
                if record.is_secondary or record.is_supplementary:
                    continue
                query_name = record.query_name
                output_number = None
                
                
                output_number = rid_to_number[query_name]
                
                output_file = dic_of_files.get(output_number, None)
                if not output_file:
                    output_file = pysam.AlignmentFile(f'{out_prefix}{output_number}.bam', "wb", template=bam)    
                    dic_of_files[output_number] = output_file               
                output_file.write(record)
        for f in dic_of_files.values():
            f.close()
    
    rid_to_number = split_ids('$ids')
    split_bam('$bef_bam', 'bef_bam_', rid_to_number)
    split_bam('$aft_bam', 'aft_bam_', rid_to_number)
    """
}

workflow splitWorkloadWorkflow {

    take:
        bef_bam
        aft_bam
        ids
        size
    main:
        splitWorkload(bef_bam, aft_bam, ids, size)
        bef_bam_split = splitWorkload.out.bef_bam.flatten()
        aft_bam_split = splitWorkload.out.aft_bam.flatten()
        ids_split = splitWorkload.out.ids.flatten()
        
    
    emit:
        bef_bam_ch = bef_bam_split
        aft_bam_ch = aft_bam_split
        ids_ch = ids_split
}
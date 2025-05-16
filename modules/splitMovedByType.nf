#!/usr/bin/env nextflow


process splitMovedByType {

    publishDir "${params.publish_dir}/split_moved_by_type", mode: 'copy'


    input:
    path moved

    output:
    path 'hap_flip.tsv', emit: hap_flip
    path 'acrocentric.tsv'
    path 'others.tsv'
    

    script:
    """
    #!/usr/bin/env python

    acrocentric_chrs = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']

    with open('$moved', 'r') as f, open('hap_flip.tsv', 'w') as hap_flip, \\
        open('acrocentric.tsv', 'w') as acrocentric, open('others.tsv', 'w') as others:
        for line in f:
            parts = line.split()
            
            bef_chr = parts[0]
            bef_chr_string_parts = bef_chr.split('_')
            bef_chrNum = bef_chr_string_parts[0]
            
            aft_chr = parts[3]
            aft_chr_string_parts = aft_chr.split('_')
            aft_chrNum = aft_chr.split('_')[0]

            if bef_chrNum in acrocentric_chrs:
                acrocentric.write(line)
                continue
            
            if len(bef_chr_string_parts) == 2 and len(aft_chr_string_parts) == 2 and \\
            bef_chrNum == aft_chrNum:
                
                bef_hap = bef_chr_string_parts[1]
                aft_hap = aft_chr_string_parts[1]
                if bef_hap != aft_hap:
                    hap_flip.write(line)
                    continue
            others.write(line)
    """
}

process combineWithTriobin {
    publishDir "${params.publish_dir}/split_moved_by_type", mode: 'copy'

    input:
    path hap_flips
    path bef_triobin, name: 'bef_triobin'
    path aft_triobin, name: 'aft_triobin'

    output:
    path 'hap_flip.triobin.tsv'
    path 'transition_counts.txt'


    script:
    """
    #!/usr/bin/env python

    from collections import defaultdict

    dic = defaultdict(list)

    def collect_into_dict(path:str, dic):
        with open(path, 'r') as f:
            for line in f:
                parts = line.split()
                readname = parts[0]
                bin = parts[1]
                pat_count = int(parts[2])
                mat_count = int(parts[3])
                pat_count1 = int(parts[4])
                mat_count1 = int(parts[5])
                #bin_specific = bin
                #if bin == 'p' or bin == 'm':
                #    if pat_count >= 21 and mat_count <=2 and bin == 'p':
                #        bin_specific = 'p'
                #    elif mat_count >= 21 and pat_count <=2 and bin == 'm':
                #        bin_specific = 'm'
                #    else:
                #        bin_specific = 'a'

                dic[readname].append((bin, pat_count, mat_count, pat_count1, mat_count1))
    
    collect_into_dict('$bef_triobin', dic)
    collect_into_dict('$aft_triobin', dic)

    with open('$hap_flips', 'r') as f, open('hap_flip.triobin.tsv', 'w') as out, open('transition_counts.txt', 'w') as transition:
        transition_counts = defaultdict(int)
        #transition_counts_specific = defaultdict(int)
        for line in f:
            readname = line.split()[6]
            bef_bin, bef_pat_count, bef_mat_count, bef_pat_count1, bef_mat_count1 = dic[readname][0]
            aft_bin, aft_pat_count, aft_mat_count, aft_pat_count1, aft_mat_count1 = dic[readname][1]
            out.write(f'{line.rstrip()}\\t{bef_pat_count}\\t{bef_mat_count}\\t{bef_pat_count1}\\t{bef_mat_count1}\\t{aft_pat_count}\\t{aft_mat_count}\\t{aft_pat_count1}\\t{aft_mat_count1}\\t{bef_bin}->{aft_bin}\\n')
            transition_counts[(bef_bin, aft_bin)] += 1
            #transition_counts_specific[(bef_bin_specific, aft_bin_specific)] += 1
        total = sum([count for count in transition_counts.values()])
        transition_counts = list(transition_counts.items())
        transition_counts.sort(reverse=True, key=lambda x: x[1])
        for (bef, aft), count in transition_counts:
            transition.write(f'{bef}->{aft}: {count}, {round(count/total*100, 2)}%\\n')
        transition.write(f'total: {total}\\n')
        #transition.write(f'higher specificity:\\n')
        #for (bef, aft), count in transition_counts_specific.items():
        #    transition.write(f'{bef}->{aft}: {count}, {round(count/total*100, 2)}%\\n')
        #transition.write(f'total: {total}\\n')
    """

}
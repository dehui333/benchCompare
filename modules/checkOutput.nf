#!/usr/bin/env nextflow


process checkOutput {
    cpus 1
    memory 5.GB

    input:
    tuple val(idx),  
        path('remaining_errors.txt'), 
        path('removed_errors.txt'), 
        path('added_errors.txt'), 
        path('beforeBench.txt'),    
        path('afterBench.txt')

    output:
    tuple val(idx), path('debug.txt')

    script:
    """
    #!/usr/bin/env python

    from collections import defaultdict

    def load_from_bench_output(path:str):
        dic = defaultdict(set)
        with open(path, 'r') as f:
            for line in f:
                parts = line.split()
                start = int(parts[1])
                end = int(parts[2])
                read_string = parts[10]
                read_string_parts = read_string.split('_') 
                readname = '_'.join(read_string_parts[:-4])
                dic[readname].add((start, end))
        return dic
    def load_single_record(path:str):
        dic = defaultdict(set)
        with open(path, 'r') as f:
            for line in f:
                parts = line.split()
                readname = parts[6]
                start = int(parts[1])
                end = int(parts[2])
                dic[readname].add((start, end))
        return dic

    def load_removed(path:str):
        dic = defaultdict(set)
        with open(path, 'r') as f:
            for line in f:
                parts = line.split()
                readname = parts[0]
                if len(parts) < 4:
                    continue
                errors_string = parts[3]
                errors_string_parts = errors_string.split('|')
                for single_error in errors_string_parts:
                    single_error_parts = single_error.split(',')
                    start = int(single_error_parts[0])
                    end = int(single_error_parts[1])
                    dic[readname].add((start, end))
        return dic
    def load_double_record(path:str):
        dic1 = defaultdict(set)
        dic2 = defaultdict(set)
        with open(path, 'r') as f:
            for line in f:
                parts = line.split()
                readname = parts[6]
                start = int(parts[1])
                end = int(parts[2])

                start2 = int(parts[8])
                end2 = int(parts[9])
                dic1[readname].add((start, end))
                dic2[readname].add((start2, end2))
        return dic1, dic2

    before_bench = load_from_bench_output('beforeBench.txt')
    after_bench = load_from_bench_output('afterBench.txt')
    added = load_single_record('added_errors.txt')
    remaining_bef, remaining_aft = load_double_record('remaining_errors.txt')
    removed = load_removed('removed_errors.txt')
    
    all_ids = set(before_bench.keys()) | set(after_bench.keys())
    with open('debug.txt', 'w') as output:
        for readname in all_ids:
            before_spans = before_bench[readname]
            after_spans = after_bench[readname]
            added_spans = added[readname]
            remaining_spans_bef = remaining_bef[readname]
            remaining_spans_aft = remaining_aft[readname]
            removed_spans = removed[readname]

            # all in remaining appear in after
            check1 = remaining_spans_aft.issubset(after_spans)
            
            # all in added appear in after
            check2 = added_spans.issubset(after_spans)
            
            # all in remaining does not appear in added and vice versa
            check3 = (len(remaining_spans_aft & added_spans) == 0)

            # num records remaining + num records added = num records afterbench
            check4 = len(remaining_spans_aft | added_spans) == len(after_spans)

            # all in removed in before 
            check5 = removed_spans.issubset(before_spans)
            
            # all remaining in beforebench
            check6 = remaining_spans_bef.issubset(before_spans)
            
            # all in removed and remaining are mutually exclusive
            check7 = (len(removed_spans & remaining_spans_bef) == 0)
            
            # num records in remaining + num records removed = num records beforebench
            check8 = len(remaining_spans_bef | removed_spans) == len(before_spans)

            # added not in before
            check9 = len(added_spans) == 0 or (not (added_spans.issubset(before_spans)))

            # removed not in after
            check10 = len(removed_spans) == 0 or (not removed_spans.issubset(after_spans))

            output.write(f'{readname}\\t{check1}\\t{check2}\\t{check3}\\t{check4}\\t{check5}\\t{check6}\\t{check7}\\t{check8}\\t{check9}\\t{check10}\\n')

    """
}
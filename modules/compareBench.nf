#!/usr/bin/env nextflow

/*
 * Take reachbench outputs before/after correction and compare them. Assuming the before and after reads
 * with the same names map to the same reference sequence.
 */


process compareBench {
    cpus 1
    memory 10.GB
    input:
    tuple val(idx), path('beforeBench.txt'), path('afterBench.txt')

    output:
    tuple val(idx), path('remaining_errors.txt'), path('removed_errors.txt'), path('added_errors.txt') 

    script:
    """
#!/usr/bin/env python

class SingleError:
    def __init__(self, ref_start:int, ref_end:int, read_pos:int, ref_seq:str, read_seq:str, is_HET:bool):
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_pos = read_pos
        self.ref_seq = ref_seq
        self.read_seq = read_seq
        self.is_HET = is_HET
        
        #self.is_HP = False
        
        #if read_seq[0] == ref_seq[0] and \
        #    (len(read_seq) >= 5 or len(ref_seq) >= 5) and \
        #    ref_seq == ref_seq[0] * len(ref_seq) and \
        #    read_seq == read_seq[0] * len(read_seq):            
        #    self.is_HP = True 
    def __repr__(self) -> str:
        return f'{self.ref_start}\\t{self.ref_end}\\t{self.read_pos}\\t{self.ref_seq}\\t{self.read_seq}\\t{"HET" if self.is_HET else "ERR"}'

    def to_string(self, ref_name:str, readname:str):
        return f'{ref_name}\\t{self.ref_start}\\t{self.ref_end}\\t{self.ref_seq}\\t{self.read_seq}\\t{"HET" if self.is_HET else "ERR"}\\t{readname}\\t{self.read_pos}'

    def compact_string(self) -> str:
        return f'{self.ref_start},{self.ref_end},{self.read_pos}'

def overlapping_pair_relationship(error1:SingleError, error2:SingleError):
    if error1.ref_start == error2.ref_start and error1.ref_end == error2.ref_end and \\
    error1.read_seq == error2.read_seq:
        return 'SAME'
    else:
        return 'DIFF'
    
class ErrorPair:
    def __init__(self, error1:SingleError, error2:SingleError):
        self.error1 = error1
        self.error2 = error2
        self.relationship = overlapping_pair_relationship(error1, error2)
    
    def __repr__(self) -> str:
        return f'{self.error1}\\t{self.error2}\\t{self.relationship}'

    def to_string(self, ref_name:str, readname:str):
        return f'{self.error1.to_string(ref_name, readname)}\\t' + \\
        f'{self.error2.ref_start}\\t{self.error2.ref_end}\\t' + \\
        f'{self.error2.ref_seq}\\t{self.error2.read_seq}\\t' + \\
        f'{"HET" if self.error2.is_HET else "ERR"}\\t{self.error2.read_pos}\\t{self.relationship}'
class ReadErrors:
    def __init__(self, read_name:str, ref_name:str) -> None:
        self.read_name = read_name
        self.ref_name = ref_name
        self.errors_list = []
    
    def add_error(self, ref_start:int, ref_end:int, read_pos:int, ref_seq:str, read_seq:str, is_HET:bool) -> None:
        self.errors_list.append(SingleError(ref_start, ref_end, read_pos, ref_seq, read_seq, is_HET))

    

    def compare(self, other):
        overlapping_pairs = []
        this_no_overlap = []
        other_no_overlap = []

        this_iter = iter(self.errors_list)
        other_iter = iter(other.errors_list)

        this_error = next(this_iter, None)
        other_error = next(other_iter, None)

        this_has_overlap = False
        other_has_overlap = False

        
        while this_error:
            # has this no other
            if not other_error:
                if this_has_overlap:
                    this_has_overlap = False
                else:
                    this_no_overlap.append(this_error)
                this_error = next(this_iter, None)
                continue
            # overlapping?
            if max(this_error.ref_start, other_error.ref_start) <= min(this_error.ref_end, other_error.ref_end):
                overlapping_pairs.append(ErrorPair(this_error, other_error))
                this_has_overlap = True
                other_has_overlap = True
                if this_error.ref_end < other_error.ref_end:
                    this_error = next(this_iter, None)
                    this_has_overlap = False
                elif this_error.ref_end > other_error.ref_end:
                    other_error = next(other_iter, None)
                    other_has_overlap = False
                else:
                    this_error = next(this_iter, None)
                    other_error = next(other_iter, None)
                    this_has_overlap = False
                    other_has_overlap = False
            else:
                if this_error.ref_end < other_error.ref_end:
                    if this_has_overlap:
                        this_has_overlap = False    
                    else:
                        this_no_overlap.append(this_error)
                    this_error = next(this_iter, None)
                    
                else:
                    if other_has_overlap:
                        other_has_overlap = False
                    else:    
                        other_no_overlap.append(other_error)
                    other_error = next(other_iter, None)
                    
        
        # has other no this
        while other_error:
            if other_has_overlap:
                other_has_overlap = False
            else:
                other_no_overlap.append(other_error)
            other_error = next(other_iter, None)
            

        #if len(overlapping_pairs) == 0:
        #    if this_error:
        #        this_no_overlap.append(this_error)
        #    if other_error:
        #        other_no_overlap.append(other_error)
        #else:
        #    if this_error and overlapping_pairs[-1][0].ref_start != this_error.ref_start:
        #        this_no_overlap.append(this_error)
        #    if other_error and overlapping_pairs[-1][1].ref_start != other_error.ref_start:
        #        other_no_overlap.append(other_error)
            

        return overlapping_pairs, this_no_overlap, other_no_overlap
        


def construct_dictionary(readbench_output:str):
    dic = {}
    with open(readbench_output, 'r') as f:
        for line in f:
            parts = line.split()
            ref_name = parts[0]
            ref_start = int(parts[1])
            ref_end = int(parts[2])
            is_HET = (parts[9] == 'HET')
            read_string = parts[10]
            read_string_parts = read_string.split('_')
            read_seq = read_string_parts[-2]
            ref_seq = read_string_parts[-3]
            read_pos = int(read_string_parts[-4])
            read_name = '_'.join(read_string_parts[:-4])
            if read_name not in dic:
                new_one = ReadErrors(read_name, ref_name)
                new_one.add_error(ref_start, ref_end, read_pos, ref_seq, read_seq, is_HET)
                dic[read_name] = new_one
            else:
                dic[read_name].add_error(ref_start, ref_end, read_pos, ref_seq, read_seq, is_HET)
    return dic

before_dict = construct_dictionary('beforeBench.txt')
after_dict = construct_dictionary('afterBench.txt')

before_set = set(before_dict.keys())
after_set = set(after_dict.keys())

intersect_set = before_set & after_set 
only_in_before = before_set - after_set
only_in_after = after_set - before_set

with open('removed_errors.txt', 'w') as cc, open('added_errors.txt', 'w') as oc, open('remaining_errors.txt', 'w') as uc:
    
    for readname in only_in_before:
        before_errors = before_dict[readname]    
        cc.write(f'{readname}\\t{before_errors.ref_name}\\t{len(before_errors.errors_list)}\\t')
        cc.write('|'.join(item.compact_string() for item in before_errors.errors_list))
        cc.write('\\n')

    for readname in only_in_after:
        after_errors = after_dict[readname]
        for item in after_errors.errors_list:
            #oc.write(f'{after_errors.ref_name}\\t{item}\\t{readname}\\n')
            oc.write(f'{item.to_string(after_errors.ref_name, readname)}\\n')
    for readname in intersect_set:
        before_errors = before_dict[readname]
        after_errors = after_dict[readname]
        
        undercorrect, correct_correct, overcorrect = before_errors.compare(after_errors)
        for item in undercorrect:
            #uc.write(f'{before_errors.ref_name}\\t{item}\\t{readname}\\n')
            uc.write(f'{item.to_string(before_errors.ref_name, readname)}\\n')
        cc.write(f'{readname}\\t{before_errors.ref_name}\\t{len(correct_correct)}\\t')
        cc.write('|'.join(item.compact_string() for item in correct_correct))
        cc.write('\\n')
        
        for item in overcorrect:
            #oc.write(f'{before_errors.ref_name}\\t{item}\\t{readname}\\n')
            oc.write(f'{item.to_string(before_errors.ref_name, readname)}\\n')
    """
}



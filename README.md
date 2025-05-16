# benchCompare

This is a Nextflow pipeline that makes use of outputs from the readbench command of [q100bench](https://github.com/nhansen/q100bench) to find remaining/added errors in a set of reads after a correction procedure. 

## Requirements

Certain tools have to be runnable via some kind of environment, e.g. local, conda, containers etc - to be configured in the config file in the repository. 

- [pysam](https://pysam.readthedocs.io/en/stable/)
- [samtools](https://www.htslib.org)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [q100bench](https://github.com/nhansen/q100bench)  

## Inputs and parameters

The inputs and parameters can be specified in the config file in the repository. Descriptions of what they are can be found in the config file.   

## Outputs

### main 
- added_errors.txt: Errors that do not exist in the uncorrected reads but exist in the corrected reads.
- remaining_errors.txt: Errors that occur at the same positions as existing errors in the uncorrected reads, or are directly adjacent. Each line contains a pair of errors before/after correction that overlap or are adjacent, 
  with the first in the uncorrected reads and the second in the corrected reads. One error overlapping with multiple will show up in multiple lines - to get the unique number of remaining errors after correction,
  find number of unique combinations of columns 7, 9 and 10 e.g. cut -f7,9,10 remaining_errors.txt | sort -u | wc -l 
  
- added_errors.txt.intersected: Added errors that occur in regions specified in the optional bed file, 'separation_bed', in the config.
- added_errors.txt.not_intersected: Added errors that do not occur in regions specified in the optional bed file, 'separation_bed', in the config.
- remaining_errors.txt.intersected: Remaining errors that occur in regions specified in the optional bed file, 'separation_bed', in the config.
- remaining_errors.txt.not_intersected: Remaining errors that do not occur in regions specified in the optional bed file, 'separation_bed', in the config.

### auxiliary 
- removed_errors.txt: Read positions from where errors were removed, details omitted.
- before_stats.txt: Contains the number of aligned bases before correction.
- after_stats.txt: Contains the number of aligned bases after correction.
- debug.txt: Results of some simple sanity check.
- moved.tsv: Information on reads whose alignments to the reference changed 'too much' after correction (>5% outside of old alignment, changed chromosome/haplotype), and are excluded from the errors analysis.
- not_moved.tsv: Information on reads whose alignments to the reference are largely overlapping before/after correction, and are included in the errors analysis.
- not_in_bam.ids: Reads that only appear in one of the before/after sets.
- unmapped_in_both.ids: Reads that are unmapped in both before/after set.
- split_moved_by_type/acrocentric.tsv: Those in not_moved.tsv that were mapped to acrocentric chromosomes before correction.
- split_moved_by_type/hap_flip.tsv: Those in not_moved.tsv that stayed in the same chromosome number but changed haplotype.
- split_moved_by_type/others.tsv: All others in not_moved.tsv that do not belong to the previous two cases.
- split_moved_by_type/hap_flip.triobin.tsv: hap_flip.tsv augmented with information from the optional yak triobin output.
- split_moved_by_type/transition_counts.txt: Counts of how often each type of triobin conclusion switch happen before/after correction. e.g. from maternal to paternal, from ambiguous to maternal    

## Brief description of steps 
1. Filter the bams if the optional bed file 'include_bed' is provided in the config.
2. Split workload into smaller chunks as specified in the config, by splitting bams etc.
3. Find reads whose alignment to the reference have low overlap between before/after correction (>5% outside of old alignment, changed chromosome/haplotype), exclude these from the error analysis.
4. Use readbench command of q100bench to find errors in reads before/after correction. Optionally use the bed file 'q100bench_bed' to limit analysis to specific regions.
5. Use the outputs from the previous step to find remaining errors - errors in corrected reads that overlap with the positions (on the reference) of errors in uncorrected reads,
   or are adjacent. The other errors are considered added errors.
6. Use the optional bed file 'separation_bed' to classify added and remaining errors in to those that fall into the regions specified in the bed file, or not.
7. Classify 'moved' reads whose alignments changed much between before/after correction into categories. Use optional yak triobin output for further analysis of cases of change in the
haplotype mapped to.   

#!/usr/bin/env nextflow

/*
 * From outputs from different splits.
 * 
 */

process mergeOutputs {
  publishDir params.publish_dir, mode: 'copy'
  cpus 1   

  input:
  path 'remaining_errors*.txt'
  path 'removed_errors*.txt'
  path 'added_errors*.txt'
  path 'moved*.tsv'
  path 'not_moved*.tsv'
  path 'unmapped_in_both*.ids'
  path 'not_in_bam*.ids'
  path 'debug*.txt'

  output:
  path 'remaining_errors.txt', emit: remaining_errors
  path 'removed_errors.txt', emit:removed_errors
  path 'added_errors.txt', emit: added_errors
  path 'moved.tsv', emit: moved
  path 'not_moved.tsv', emit: not_moved
  path 'unmapped_in_both.ids', emit: unmapped
  path 'not_in_bam.ids', emit: not_in_bam
  path 'debug.txt', emit: debug
  script:
  """
  cat remaining_errors*.txt > remaining_errors.txt
  cat removed_errors*.txt > removed_errors.txt
  cat added_errors*.txt > added_errors.txt
  cat moved*.tsv > moved.tsv
  cat not_moved*.tsv > not_moved.tsv
  cat unmapped_in_both*.ids > unmapped_in_both.ids
  cat not_in_bam*.ids > not_in_bam.ids
  cat debug*.txt > debug.txt
  """
}

process getStats {
  publishDir params.publish_dir, mode: 'copy'
  cpus 1

  input:
  path 'stats*.txt'
  val output_prefix

  output:
  path "${output_prefix}_stats.txt"

  script:
  """
  #!/usr/bin/env python
  import glob

  def extract_aligned_bases_num(path:str) -> int:
    with open(path, 'r') as f:
      return int(f.readline().split()[3]) 

  input_stats = glob.glob('stats*')
  with open('${output_prefix}_stats.txt', 'w') as out:
    out.write(f'Total number of aligned bases: {sum(map(extract_aligned_bases_num, input_stats))}\\n')
  """
}

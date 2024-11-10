# BIOINF_TOOLS

Bioinf_tools contains tool for filtering FASTQ sequences, tool for making simple transformations of RNA and DNA 
sequences, tool for converting multiline FASTA sequences to one-line FASTA sequences, and tool for extracting names of
top matching protein for each query sequence from a BLAST result file.

Authors: Kozulina Maria

## Running instructions

To run bioinf_tools use the following syntax:

```
run_dna_rna_tools(*args)
```
```
filter_fastq(input_fastq: str | Path,
             output_fastq: str | Path,
             gc_bounds: tuple[Real, Real] | list[Real, Real] | Real = (0, 100),
             length_bounds: tuple[int, int] | list[int, int] | int = (0, 2 ** 32),
             quality_threshold: Real = 0,
             rewrite: bool = False)
```
To run bio_files_processor use the following syntax:
```
convert_multiline_fasta_to_oneline(input_fasta: str | Path, output_fasta: str | Path | None = None)
```
```
parse_blast_output(input_file: str | Path, output_file: str | Path)
```

## Features
### filter_fastq
The program uses file with FASTQ data, path for file with results, numbers for length range, GC-content and quality-threshold for input.  
The program filters FASTQ sequences based on GC content, length of the sequence, and mean quality of the reads (PHRED33).  The length of the sequence has default range(0, 2^32). The GC content has default range (0, 100) in percentage. 
Default value for quality threshold is 0.  The program returns writes a new file with FASTQ data for the sequences that meet all the criteria from input.
### run_rna_dna_tools
The program uses 1 sequence or number of sequences of RNA and DNA and applies a given function to the sequence/sequences.
The supported functions are:
* function "transcribe" returns transcribed sequence
* function "reverse" returns reversed sequence
* function "complement" returns complemented sequence
* function "reverse_complement" returns complemented and reversed sequence
* function "find_stop" returns the full sequence if no stop codon is present. If a stop codon is found, it returns the sequence up to and including the first stop codon. 
* function "find_start"  returns the sequence starting from the first start codon. If no start codon is found, it returns an empty string.
* function "check_seq" checks if the string from input is DNA or RNA sequence.
### convert_multiline_fasta_to_oneline
The program uses file with FASTA data and creates a file with all the sequences from the data written in one line.
### parse_blast_output
The program extracts the description of the best matching sequence from the Description section for each query.  The results are saved in a new text file in one column, sorted alphabetically.

## Examples
```python
run_dna_rna_tools('AtAcGGC','ccAUcTc' 'transcribe') -> ['AuAcGGC', 'ccATctc']
```
```python
run_dna_rna_tools('AUcgCGCCauaAccGccuUaa','find_stop') -> 'AUcgCGCCauaA'
```
```python
filter_fastq('../input_file', '../output_file', gc_bounds = (20,50), length_bounds = (10,28), 
            quality_threshold = 30, rewrite = True)
```

## Project structure
```
-/
 |- README.md
 |- bioinf_tools.py
 |- bio_files_processor.py
 |- config.py
 |- src/
       |- dna_rna_tools.py - utils for run_dna_rna_tools.
       |- filter_fastq_tools.py - utils for filter_fastq.
```

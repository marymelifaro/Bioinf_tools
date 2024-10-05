# BIOINF_TOOLS

Bioinf_tools is a tool for filtering FASTQ sequences and making simple transformations of RNA and DNA sequences. 

Authors: Kozulina Maria

## Running instructions

To run bioinf_tools use the following syntax:

```
run_dna_rna_tools(*args)
```
```
filter_fastq((seqs: dict[str, tuple[[str, str]]],
                 gc_bounds: tuple[number, number] | number = (0, 100),
                 length_bounds: tuple[int, int] | int = (0, 2 ** 32),
                 quality_threshold: number = 0) -> dict[str, tuple[[str, str]]]:)
```
## Features
### filter_fastq
The program uses dictionary with FASTQ data (name of sequence as key, list with sequence and its quality as value), numbers for length range, GC-content and quality-threshold for input.  
The program filters FASTQ sequences based on GC content, length of the sequence, and mean quality of the reads (PHRED33).
The length of the sequence has default range(0, 2^32). The GC content has default range (0, 100) in percentage. Default value for quality threshold is 0.  The program returns new dictionary with FASTQ data for the sequences that meet all the criteria from input.
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

## Examples
```python
run_dna_rna_tools('AtAcGGC','ccAUcTc' 'transcribe') -> ['AuAcGGC', 'ccATctc']
```
```python
run_dna_rna_tools('AUcgCGCCauaAccGccuUaa','find_stop') -> 'AUcgCGCCauaA'
```
```python
filter_fastq({'@SRX079804': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG', 
                             'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'), 
              '@SRX079809': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT', 
                             'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE')},(20,50),(10,80), 30)) -> 
{'@SRX079809': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT', 'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE')}
```

## Project structure
```
-/
 |- README.md
 |- bioinf_tools.py - main script
 |- src/
       |- dna_rna_tools.py - utils for run_dna_rna_tools.
       |- filter_fastq_tools.py - utils for filter_fastq.
```

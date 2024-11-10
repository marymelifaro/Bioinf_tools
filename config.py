OUTPUT_FASTA_MULTILINE = '_output'

STOP_CODONS = ['UAA', 'UAG', 'UGA']

DICT_TRANS = {'A': 'A', 'a': 'a', 'T': 'U', 't': 'u', 'U': 'T',
              'u': 't', 'G': 'G', 'g': 'g', 'C': 'C', 'c': 'c'}

DICT_COMPL = {'A': 'T', 'a': 't', 'T': 'A', 't': 'a', 'U': 'A',
              'u': 'a', 'G': 'C', 'g': 'c', 'C': 'G', 'c': 'g'}

START_CODON = 'AUG'
DNA_NUCLEOTIDES = set('ATCG')
RNA_NUCLEOTIDES = set('AUCG')
FORBIDDEN_NUCLEOTIDES = set('TU')

PHRED33_SHIFT = 33

FASTQ_SAMPLE_ROWS = 4

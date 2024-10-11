from typing import Any

from config import (DICT_COMPL, DICT_TRANS, DNA_NUCLEOTIDES, FORBIDDEN_NUCLEOTIDES, RNA_NUCLEOTIDES, START_CODON,
                    STOP_CODONS)


def transcribe(seq: str) -> str:
    return ''.join([DICT_TRANS[x] for x in seq])


def reverse(seq: str) -> str:
    return seq[::-1]


def complement(seq: str) -> str:
    return ''.join([DICT_COMPL[x] for x in seq])


def reverse_complement(seq: str) -> str:
    return reverse(complement(seq))


def find_stop(seq: str) -> str:
    for j in range(0, len(seq), 3):
        if seq[j:j+3].upper() in STOP_CODONS:
            return seq[:j+3]
    return seq


def find_start(seq: str) -> str:
    for j in range(0, len(seq), 3):
        if seq[j:j+3].upper() == START_CODON:
            return seq[j:]
    return ''


def check_seq(seq: Any):
    if not isinstance(seq, str):
        raise TypeError(f"sequence must be str type. Now type(seq) = {type(seq)}")
    seq = seq.upper()
    unique_nucleotides = set(seq)
    if len(unique_nucleotides) > 4:
        raise ValueError("sequence must have not more than 4 unique symbols. "
                         f"Now unique symbols are: {unique_nucleotides}")

    if FORBIDDEN_NUCLEOTIDES.issubset(unique_nucleotides):
        raise ValueError("sequence must consist of DNA nucleotides or RNA nucleotides. "
                         "Now the sequence contains both 'T' and 'U'. ")

    if not unique_nucleotides.issubset(DNA_NUCLEOTIDES) and not unique_nucleotides.issubset(RNA_NUCLEOTIDES):
        raise ValueError(f"sequence must consist of DNA nucleotides or RNA nucleotides. "
                         f"Now unique nucleotides are {unique_nucleotides}")


DNA_RNA_FUNCTIONS = {
    'transcribe': transcribe,
    'reverse': reverse,
    'complement': complement,
    'reverse_complement': reverse_complement,
    'find_stop': find_stop,
    'find_start': find_start,
}

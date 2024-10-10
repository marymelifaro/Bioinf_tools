from numbers import Real
from src.dna_rna_tools import DNA_RNA_FUNCTIONS, check_seq
from src.filter_fastq_tools import filter_gc, filter_quality, filter_length


def run_dna_rna_tools(*args) -> str | list[str]:
    """
    Applies a given function to one or more nucleotide sequences.

    :param seqs: A string or a list of strings representing nucleotide sequences.
    :param function_name: The function to apply to the nucleotide sequence(s) (e.g., 'transcribe', 'complement').
    :return: A string or a list of strings with nucleotide sequences after applying the specified function.
    """
    *seqs, function_name = args
    res = []
    for seq in seqs:
        check_seq(seq)
        if function_name in DNA_RNA_FUNCTIONS:
            res.append(DNA_RNA_FUNCTIONS[function_name](seq))
        else:
            raise ValueError(f"function {function_name} is not supported")
    if len(res) == 1:
        return res[0]
    return res


def filter_fastq(seqs: dict[str, tuple[[str, str]]],
                 gc_bounds: tuple[Real, Real] | list[Real, Real] | Real = (0, 100),
                 length_bounds: tuple[int, int] | list[int, int] | int = (0, 2 ** 32),
                 quality_threshold: Real = 0) -> dict[str, tuple[[str, str]]]:
    """
    Filters FASTQ sequences based on GC content, length, and quality.

    :param seqs: Dictionary with FASTQ sequences. The key is a string representing the name of the sequence.
                 The value is a tuple consisting of the sequence itself (as a string) and
                 its quality scores (as a string).
    :param gc_bounds: Range of values for GC content percentage in the sequence. A tuple with two numbers defines
                      the range, while a single number represents the upper limit. Default is (0, 100).
    :param length_bounds: Range for the length of the sequence. A tuple with two numbers represents a range.
                          Default value is (0, 2^32).
    :param quality_threshold: Threshold value for the average quality of the sequence. Reads with an average quality
                              below this threshold will be filtered out. Default value is 0.
    :return: New dictionary with FASTQ sequences that meet all the filtering criteria
            (quality threshold, length bounds, GC bounds).
    """
    if isinstance(gc_bounds, list):
        gc_bounds = tuple(gc_bounds)
    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    if isinstance(length_bounds, list):
        length_bounds = tuple(length_bounds)
    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)

    filtered_seqs = {}
    for name, (seq, quality) in seqs.items():
        check_seq(seq)
        if not (filter_length(seq, length_bounds) and filter_gc(seq, gc_bounds)
                and filter_quality(quality, quality_threshold)):
            continue
        filtered_seqs[name] = seqs[name]
    return filtered_seqs

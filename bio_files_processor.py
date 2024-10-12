from pathlib import Path

from config import OUTPUT_FASTA_MULTILINE


def convert_multiline_fasta_to_oneline(input_fasta: str | Path, output_fasta: str | Path | None = None):
    """
    Converts multiline FASTA sequence to 1 line.

    :param input_fasta: string with the path to file with FASTA sequences.
    :param output_fasta: string with the path to file with FASTA sequences, written in 1 line.
                        If output_fasta path is None, function use input_fasta path as template for output_fasta.
    """


    input_fasta = Path(input_fasta)
    if output_fasta is None:
        output_fasta = input_fasta.parent / (input_fasta.stem + OUTPUT_FASTA_MULTILINE + input_fasta.suffix)

    Path(output_fasta).parent.mkdir(exist_ok=True)

    is_first_line = True
    with open(input_fasta, 'r') as in_file, open(output_fasta, 'w') as out_file:
        for line in in_file:
            if not line.startswith('>'):
                line = line.strip()
                out_file.write(line)
            else:
                if not is_first_line:
                    out_file.write('\n')
                else:
                    is_first_line = False
                out_file.write(line)


def parse_blast_output(input_file: str | Path, output_file: str | Path):
    """
    Reads a BLAST results text file, extracts the description of the best match (first hit) for each query sequence
    and writes the extracted protein names into the output file. Protein names are sorted alphabetically and are
    written in one column.

    :param input_file: string with the path to file with blast results.
    :param output_file: string with the path to output_file.
                        If output_file path is None, function use input_file path as template for output_file.
    """
    result = []
    with open(input_file, 'r') as in_file:
        cnt = -2 ** 32
        for line in in_file:
            if line.startswith('Query #'):
                cnt = 0
            cnt += 1
            if cnt == 5:
                descr_length = line.find('Name')
            if cnt == 6:
                result.append(line[:descr_length].strip())
    result = sorted(result)

    Path(output_file).parent.mkdir(exist_ok=True)
    with open(output_file, 'w') as file:
        for item in result:
            file.write(item)
            file.write('\n')

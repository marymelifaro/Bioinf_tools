from pathlib import Path

from config import OUTPUT_FASTA_MULTILINE


def convert_multiline_fasta_to_oneline(input_fasta: str | Path, output_fasta: str | Path | None = None):
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

import argparse
import math
import os
import pandas as pd
from Bio import SeqIO


NUCLEOTIDE_TO_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
CSV_COLUMNS = ["A", "C", "G", "T"]


def positional_encode(sequence: str, alpha: float = 0.05, period: float = 10.0):
    length = len(sequence)
    if length == 0:
        return []

    cls_token = [
        sequence.count("A") / length,
        sequence.count("C") / length,
        sequence.count("G") / length,
        sequence.count("T") / length,
    ]

    encoded_matrix = [cls_token]

    for position, nucleotide in enumerate(sequence, start=1):
        if nucleotide not in NUCLEOTIDE_TO_INDEX:
            continue

        one_hot = [0.0, 0.0, 0.0, 0.0]
        one_hot[NUCLEOTIDE_TO_INDEX[nucleotide]] = 1.0

        angle = 2 * math.pi * position / period
        sin_val = math.sin(angle)
        cos_val = math.cos(angle)
        pos_enc = [sin_val, cos_val, -sin_val, cos_val]

        encoded_matrix.append(
            [value + alpha * offset for value, offset in zip(one_hot, pos_enc)]
        )

    return encoded_matrix


def sanitize_filename(name: str) -> str:
    sanitized = "".join(c for c in name if c.isalnum() or c in {"_", "-"}).rstrip()
    return sanitized or "sequence"


def process_fasta_file(input_file: str, output_base_dir: str, alpha: float, period: float):
    file_name = os.path.splitext(os.path.basename(input_file))[0]
    output_dir = os.path.join(output_base_dir, file_name)
    os.makedirs(output_dir, exist_ok=True)

    processed_count = 0
    total_count = 0

    for record in SeqIO.parse(input_file, "fasta"):
        total_count += 1
        name = record.id
        sequence = str(record.seq).upper()

        encoded = positional_encode(sequence, alpha=alpha, period=period)
        if not encoded:
            continue

        output_file = os.path.join(output_dir, f"{sanitize_filename(name)}.csv")
        pd.DataFrame(encoded, columns=CSV_COLUMNS).to_csv(output_file, index=False)
        processed_count += 1

    return processed_count, total_count


def process_all_files(input_dir: str, output_dir: str, alpha: float, period: float):
    os.makedirs(output_dir, exist_ok=True)

    fasta_files = sorted(
        file_name for file_name in os.listdir(input_dir) if file_name.endswith(".txt")
    )

    if not fasta_files:
        print(f"No .txt files found in {input_dir}")
        return

    print(f"Found {len(fasta_files)} files")

    total_processed = 0
    total_sequences = 0

    for index, file_name in enumerate(fasta_files, start=1):
        input_path = os.path.join(input_dir, file_name)
        print(f"[{index}/{len(fasta_files)}] Processing {file_name}")

        try:
            processed_count, total_count = process_fasta_file(
                input_path,
                output_dir,
                alpha=alpha,
                period=period,
            )
            total_processed += processed_count
            total_sequences += total_count
            print(f"  sequences: {processed_count}/{total_count}")
        except Exception as error:
            print(f"  failed: {error}")

    print(f"Files processed: {len(fasta_files)}")
    print(f"Sequences processed: {total_processed}/{total_sequences}")



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--period", type=float, default=10.0)
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.isdir(args.input_dir):
        raise FileNotFoundError(f"Input directory does not exist: {args.input_dir}")

    process_all_files(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        alpha=args.alpha,
        period=args.period,
    )


if __name__ == "__main__":
    main()

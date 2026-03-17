"""Utility functions for sequence manipulation and validation."""

import re
from typing import Optional

VALID_DNA_BASES = set("ACGTacgt")
COMPLEMENT_MAP = str.maketrans("ACGTacgt", "TGCAtgca")


def validate_sequence(sequence: str) -> bool:
    """Validate that a string contains only valid DNA bases (A, C, G, T).

    Args:
        sequence: DNA sequence string.

    Returns:
        True if valid, False otherwise.
    """
    if not sequence:
        return False
    return all(c in VALID_DNA_BASES for c in sequence)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        sequence: DNA sequence string (A, C, G, T).

    Returns:
        Reverse complement string.

    Raises:
        ValueError: If sequence contains invalid bases.
    """
    if not validate_sequence(sequence):
        raise ValueError(f"Invalid DNA sequence: {sequence}")
    return sequence.translate(COMPLEMENT_MAP)[::-1]


def gc_content(sequence: str) -> float:
    """Calculate GC content as a fraction (0.0 to 1.0).

    Args:
        sequence: DNA sequence string.

    Returns:
        GC content fraction.

    Raises:
        ValueError: If sequence is empty or invalid.
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")
    seq_upper = sequence.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    return gc_count / len(seq_upper)


def melt_temp_basic(sequence: str) -> float:
    """Calculate melting temperature using the Wallace rule.

    For short oligonucleotides (<14 bp): Tm = 2(A+T) + 4(G+C)
    For longer sequences: Tm = 64.9 + 41*(G+C - 16.4) / N

    Args:
        sequence: DNA sequence string.

    Returns:
        Estimated melting temperature in °C.
    """
    seq_upper = sequence.upper()
    n = len(seq_upper)
    gc = seq_upper.count("G") + seq_upper.count("C")
    at = seq_upper.count("A") + seq_upper.count("T")

    if n < 14:
        return float(2 * at + 4 * gc)
    else:
        return 64.9 + 41.0 * (gc - 16.4) / n


def self_complementarity_score(sequence: str) -> int:
    """Calculate a self-complementarity score for a sequence.

    Checks all possible alignments of the sequence against its own
    reverse complement to find the maximum number of contiguous
    base-pair matches.

    Args:
        sequence: DNA sequence string.

    Returns:
        Maximum contiguous self-complementary base-pair count.
    """
    seq_upper = sequence.upper()
    rc = reverse_complement(seq_upper)
    n = len(seq_upper)
    max_run = 0

    for offset in range(-(n - 1), n):
        run = 0
        for i in range(n):
            j = i - offset
            if 0 <= j < n and seq_upper[i] == rc[j]:
                run += 1
                max_run = max(max_run, run)
            else:
                run = 0
    return max_run


def has_homopolymer_run(sequence: str, max_run: int = 4) -> bool:
    """Check if the sequence has a homopolymer run exceeding max_run.

    Args:
        sequence: DNA sequence string.
        max_run: Maximum allowed homopolymer length.

    Returns:
        True if a homopolymer run exceeds max_run.
    """
    seq_upper = sequence.upper()
    for base in "ACGT":
        pattern = base * (max_run + 1)
        if pattern in seq_upper:
            return True
    return False

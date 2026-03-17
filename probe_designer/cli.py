"""Command-line interface for the Hybrid Capture Probe Designer."""

import argparse
import sys
from typing import List, Optional

from probe_designer.core.designer import ProbeDesigner


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (defaults to sys.argv[1:]).

    Returns:
        Parsed namespace.
    """
    parser = argparse.ArgumentParser(
        prog="probe-designer",
        description=(
            "Hybrid Capture Probe Designer - Design ultra-short capture probes "
            "for DNA targets (<40 bp) with optional LNA incorporation."
        ),
    )

    parser.add_argument(
        "target",
        help="Target DNA sequence (A/C/G/T only) or path to a FASTA file.",
    )
    parser.add_argument(
        "-l", "--probe-length",
        type=int,
        default=None,
        help="Fixed probe length in bp. If not set, tries multiple lengths.",
    )
    parser.add_argument(
        "-n", "--top-n",
        type=int,
        default=10,
        help="Number of top probes to return (default: 10).",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=15,
        help="Minimum probe length in bp (default: 15).",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=40,
        help="Maximum probe length in bp (default: 40).",
    )
    parser.add_argument(
        "--min-gc",
        type=float,
        default=0.30,
        help="Minimum GC content fraction (default: 0.30).",
    )
    parser.add_argument(
        "--max-gc",
        type=float,
        default=0.70,
        help="Maximum GC content fraction (default: 0.70).",
    )
    parser.add_argument(
        "--min-tm",
        type=float,
        default=45.0,
        help="Minimum melting temperature in °C (default: 45.0).",
    )
    parser.add_argument(
        "--max-tm",
        type=float,
        default=75.0,
        help="Maximum melting temperature in °C (default: 75.0).",
    )
    parser.add_argument(
        "--optimal-tm",
        type=float,
        default=60.0,
        help="Target melting temperature in °C (default: 60.0).",
    )
    parser.add_argument(
        "--na-conc",
        type=float,
        default=0.05,
        help="Sodium concentration in M (default: 0.05).",
    )
    parser.add_argument(
        "--mg-conc",
        type=float,
        default=0.0,
        help="Magnesium concentration in M (default: 0.0).",
    )
    parser.add_argument(
        "--no-lna",
        action="store_true",
        help="Disable LNA incorporation.",
    )
    parser.add_argument(
        "--single-strand",
        action="store_true",
        help="Design probes on sense strand only.",
    )
    parser.add_argument(
        "--tiling-step",
        type=int,
        default=1,
        help="Tiling step size in bp (default: 1).",
    )

    return parser.parse_args(argv)


def _read_fasta(filepath: str) -> str:
    """Read the first sequence from a FASTA file.

    Args:
        filepath: Path to FASTA file.

    Returns:
        DNA sequence string.
    """
    sequence_parts = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence_parts:
                    break  # only read first sequence
                continue
            sequence_parts.append(line)
    return "".join(sequence_parts)


def main(argv: Optional[List[str]] = None) -> int:
    """Main entry point for the CLI.

    Args:
        argv: Argument list (defaults to sys.argv[1:]).

    Returns:
        Exit code (0 for success).
    """
    args = parse_args(argv)

    # Determine if target is a file or sequence
    target = args.target
    if target.endswith((".fa", ".fasta", ".fna")):
        try:
            target = _read_fasta(target)
        except FileNotFoundError:
            print(f"Error: File not found: {args.target}", file=sys.stderr)
            return 1
        except (IOError, OSError) as e:
            print(f"Error reading file: {e}", file=sys.stderr)
            return 1

    # Create designer
    designer = ProbeDesigner(
        min_probe_length=args.min_length,
        max_probe_length=args.max_length,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        min_tm=args.min_tm,
        max_tm=args.max_tm,
        optimal_tm=args.optimal_tm,
        use_lna=not args.no_lna,
        na_conc=args.na_conc,
        mg_conc=args.mg_conc,
        both_strands=not args.single_strand,
        tiling_step=args.tiling_step,
    )

    try:
        summary = designer.design_probes_summary(
            target,
            probe_length=args.probe_length,
            top_n=args.top_n,
        )
        print(summary)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

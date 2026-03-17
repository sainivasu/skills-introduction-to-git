"""Main probe design engine for hybrid capture of ultra-short DNA targets.

Generates candidate capture probes from target sequences, evaluates them
using thermodynamic criteria, and optionally incorporates LNA modifications
to enhance binding to short targets (<40 bp).

Design strategy:
  1. Generate tiling probes across the target region
  2. Filter by GC content, homopolymer runs, and self-complementarity
  3. Calculate thermodynamic properties (Tm, ΔG)
  4. Optionally incorporate LNA to boost Tm for ultra-short probes
  5. Rank and return the best probes
"""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from probe_designer.core.thermodynamics import ThermodynamicCalculator
from probe_designer.core.lna import LNADesigner
from probe_designer.core.utils import (
    gc_content,
    has_homopolymer_run,
    reverse_complement,
    self_complementarity_score,
    validate_sequence,
)


@dataclass
class ProbeCandidate:
    """A candidate capture probe with its properties.

    Attributes:
        sequence: Probe DNA sequence (5' to 3').
        start: 0-based start position on the target.
        end: 0-based end position (exclusive) on the target.
        strand: '+' for sense, '-' for antisense.
        gc: GC content fraction.
        tm: Predicted melting temperature in °C.
        dg_37: Free energy at 37°C in kcal/mol.
        dh: Enthalpy in kcal/mol.
        ds: Entropy in cal/(mol·K).
        self_comp: Self-complementarity score.
        lna_positions: LNA-modified positions (0-based).
        lna_sequence: Sequence with LNA notation.
        score: Overall quality score (higher = better).
        warnings: List of any quality warnings.
    """

    sequence: str
    start: int
    end: int
    strand: str = "+"
    gc: float = 0.0
    tm: float = 0.0
    dg_37: float = 0.0
    dh: float = 0.0
    ds: float = 0.0
    self_comp: int = 0
    lna_positions: List[int] = field(default_factory=list)
    lna_sequence: str = ""
    score: float = 0.0
    warnings: List[str] = field(default_factory=list)


class ProbeDesigner:
    """Design hybrid capture probes for ultra-short DNA targets.

    Args:
        min_probe_length: Minimum probe length in bp (default 15).
        max_probe_length: Maximum probe length in bp (default 40).
        min_gc: Minimum GC content fraction (default 0.30).
        max_gc: Maximum GC content fraction (default 0.70).
        min_tm: Minimum melting temperature in °C (default 45.0).
        max_tm: Maximum melting temperature in °C (default 75.0).
        optimal_tm: Target melting temperature in °C (default 60.0).
        max_self_comp: Maximum self-complementarity score (default 6).
        max_homopolymer: Maximum homopolymer run length (default 4).
        use_lna: Whether to incorporate LNA modifications (default True).
        na_conc: Sodium concentration in M (default 0.05).
        mg_conc: Magnesium concentration in M (default 0.0).
        oligo_conc: Probe concentration in M (default 250e-9).
        both_strands: Design probes on both strands (default True).
        tiling_step: Step size for tiling probes (default 1).
    """

    def __init__(
        self,
        min_probe_length: int = 15,
        max_probe_length: int = 40,
        min_gc: float = 0.30,
        max_gc: float = 0.70,
        min_tm: float = 45.0,
        max_tm: float = 75.0,
        optimal_tm: float = 60.0,
        max_self_comp: int = 6,
        max_homopolymer: int = 4,
        use_lna: bool = True,
        na_conc: float = 0.05,
        mg_conc: float = 0.0,
        oligo_conc: float = 250e-9,
        both_strands: bool = True,
        tiling_step: int = 1,
    ) -> None:
        self.min_probe_length = min_probe_length
        self.max_probe_length = max_probe_length
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.optimal_tm = optimal_tm
        self.max_self_comp = max_self_comp
        self.max_homopolymer = max_homopolymer
        self.use_lna = use_lna
        self.both_strands = both_strands
        self.tiling_step = tiling_step

        self.thermo = ThermodynamicCalculator(
            na_conc=na_conc,
            mg_conc=mg_conc,
            oligo_conc=oligo_conc,
        )
        self.lna_designer = LNADesigner()

    def _generate_tiles(
        self, target: str, probe_length: int
    ) -> List[Tuple[str, int, int, str]]:
        """Generate tiling probe sequences from a target.

        Args:
            target: Target DNA sequence.
            probe_length: Length of each tile.

        Returns:
            List of (sequence, start, end, strand) tuples.
        """
        tiles = []
        n = len(target)

        if probe_length > n:
            # If probe is longer than target, use the full target
            tiles.append((target, 0, n, "+"))
            if self.both_strands:
                rc = reverse_complement(target)
                tiles.append((rc, 0, n, "-"))
            return tiles

        for start in range(0, n - probe_length + 1, self.tiling_step):
            end = start + probe_length
            seq = target[start:end]
            tiles.append((seq, start, end, "+"))

            if self.both_strands:
                rc = reverse_complement(seq)
                tiles.append((rc, start, end, "-"))

        return tiles

    def _evaluate_probe(
        self, sequence: str, start: int, end: int, strand: str
    ) -> Optional[ProbeCandidate]:
        """Evaluate a single probe candidate.

        Args:
            sequence: Probe DNA sequence.
            start: Start position on target.
            end: End position on target.
            strand: '+' or '-'.

        Returns:
            ProbeCandidate if it passes filters, None otherwise.
        """
        warnings: List[str] = []
        seq_upper = sequence.upper()

        # Basic validation
        if not validate_sequence(seq_upper):
            return None

        # GC content filter
        gc = gc_content(seq_upper)
        if gc < self.min_gc:
            return None
        if gc > self.max_gc:
            return None

        # Homopolymer filter
        if has_homopolymer_run(seq_upper, self.max_homopolymer):
            return None

        # Self-complementarity
        self_comp = self_complementarity_score(seq_upper)
        if self_comp > self.max_self_comp:
            return None

        # Thermodynamic calculation (without LNA first)
        thermo = self.thermo.calculate(seq_upper)
        tm = thermo["Tm"]
        dg_37 = thermo["dG_37"]
        dh = thermo["dH"]
        ds = thermo["dS"]

        lna_positions: List[int] = []
        lna_sequence = seq_upper

        # Try LNA incorporation if Tm is below target and LNA is enabled
        if self.use_lna and tm < self.optimal_tm:
            lna_positions = self.lna_designer.design_lna_pattern(seq_upper)
            thermo_lna = self.thermo.calculate(seq_upper, lna_positions)
            tm = thermo_lna["Tm"]
            dg_37 = thermo_lna["dG_37"]
            dh = thermo_lna["dH"]
            ds = thermo_lna["dS"]
            lna_sequence = self.lna_designer.format_lna_sequence(
                seq_upper, lna_positions
            )

            validation = self.lna_designer.validate_lna_pattern(
                seq_upper, lna_positions
            )
            warnings.extend(validation["warnings"])

        # Check Tm bounds after potential LNA incorporation
        if tm < self.min_tm:
            warnings.append(f"Tm ({tm:.1f}°C) below minimum ({self.min_tm}°C)")
        if tm > self.max_tm:
            warnings.append(f"Tm ({tm:.1f}°C) above maximum ({self.max_tm}°C)")

        # Score the probe
        score = self._score_probe(gc, tm, dg_37, self_comp, len(lna_positions))

        return ProbeCandidate(
            sequence=seq_upper,
            start=start,
            end=end,
            strand=strand,
            gc=round(gc, 3),
            tm=tm,
            dg_37=dg_37,
            dh=dh,
            ds=ds,
            self_comp=self_comp,
            lna_positions=lna_positions,
            lna_sequence=lna_sequence,
            score=round(score, 3),
            warnings=warnings,
        )

    def _score_probe(
        self,
        gc: float,
        tm: float,
        dg_37: float,
        self_comp: int,
        num_lna: int,
    ) -> float:
        """Calculate an overall quality score for a probe.

        Score components:
          - Tm proximity to optimal (weighted heavily)
          - GC content near 0.50 (ideal)
          - Favorable ΔG (more negative = more stable)
          - Low self-complementarity
          - Penalty for excessive LNA usage

        Args:
            gc: GC content fraction.
            tm: Melting temperature in °C.
            dg_37: Free energy at 37°C.
            self_comp: Self-complementarity score.
            num_lna: Number of LNA bases.

        Returns:
            Quality score (higher = better).
        """
        score = 0.0

        # Tm component: penalize deviation from optimal
        tm_diff = abs(tm - self.optimal_tm)
        score += max(0, 40.0 - tm_diff * 2.0)

        # GC component: prefer 40-60%
        gc_diff = abs(gc - 0.50)
        score += max(0, 20.0 - gc_diff * 40.0)

        # Stability component: more negative ΔG is better (up to a point)
        # Typical range: -5 to -25 kcal/mol for short probes
        if dg_37 < 0:
            score += min(20.0, abs(dg_37))

        # Self-complementarity penalty
        score -= self_comp * 1.5

        # Moderate LNA penalty (prefer fewer modifications when possible)
        score -= num_lna * 0.5

        return max(0.0, score)

    def design_probes(
        self,
        target_sequence: str,
        probe_length: Optional[int] = None,
        top_n: int = 10,
    ) -> List[ProbeCandidate]:
        """Design capture probes for an ultra-short DNA target.

        Main entry point for probe design. Generates candidate probes,
        evaluates them, and returns the top-ranked candidates.

        Args:
            target_sequence: Target DNA sequence to capture.
            probe_length: Fixed probe length. If None, tries multiple
                lengths within the allowed range.
            top_n: Number of top probes to return.

        Returns:
            List of ProbeCandidate objects, sorted by score (best first).

        Raises:
            ValueError: If target sequence is invalid.
        """
        target = target_sequence.upper().strip()

        if not validate_sequence(target):
            raise ValueError(
                f"Invalid target sequence. Only A, C, G, T bases are allowed."
            )

        if len(target) < self.min_probe_length:
            raise ValueError(
                f"Target sequence ({len(target)} bp) is shorter than "
                f"minimum probe length ({self.min_probe_length} bp)."
            )

        candidates: List[ProbeCandidate] = []

        if probe_length is not None:
            # Fixed probe length
            lengths = [probe_length]
        else:
            # Try multiple lengths, capped at target length
            max_len = min(self.max_probe_length, len(target))
            min_len = min(self.min_probe_length, max_len)
            lengths = list(range(min_len, max_len + 1))

        for pl in lengths:
            tiles = self._generate_tiles(target, pl)
            for seq, start, end, strand in tiles:
                probe = self._evaluate_probe(seq, start, end, strand)
                if probe is not None:
                    candidates.append(probe)

        # Sort by score descending
        candidates.sort(key=lambda p: -p.score)

        # Return top N
        return candidates[:top_n]

    def design_probes_summary(
        self,
        target_sequence: str,
        probe_length: Optional[int] = None,
        top_n: int = 10,
    ) -> str:
        """Design probes and return a formatted summary string.

        Convenience method that wraps design_probes() with formatted output.

        Args:
            target_sequence: Target DNA sequence.
            probe_length: Optional fixed probe length.
            top_n: Number of top probes.

        Returns:
            Formatted summary string.
        """
        probes = self.design_probes(target_sequence, probe_length, top_n)

        lines = []
        lines.append("=" * 72)
        lines.append("HYBRID CAPTURE PROBE DESIGN RESULTS")
        lines.append("=" * 72)
        lines.append(f"Target: {target_sequence[:50]}{'...' if len(target_sequence) > 50 else ''}")
        lines.append(f"Target length: {len(target_sequence)} bp")
        lines.append(f"Probes designed: {len(probes)}")
        lines.append(f"LNA enabled: {self.use_lna}")
        lines.append("")

        if not probes:
            lines.append("No valid probes found. Try adjusting parameters.")
            return "\n".join(lines)

        for i, p in enumerate(probes, 1):
            lines.append(f"--- Probe #{i} ---")
            lines.append(f"  Sequence:    {p.sequence}")
            if p.lna_positions:
                lines.append(f"  LNA format:  {p.lna_sequence}")
                lines.append(f"  LNA count:   {len(p.lna_positions)}")
            lines.append(f"  Position:    {p.start}-{p.end} ({p.strand} strand)")
            lines.append(f"  Length:      {len(p.sequence)} bp")
            lines.append(f"  GC content:  {p.gc * 100:.1f}%")
            lines.append(f"  Tm:          {p.tm:.1f}°C")
            lines.append(f"  ΔG (37°C):   {p.dg_37:.2f} kcal/mol")
            lines.append(f"  ΔH:          {p.dh:.2f} kcal/mol")
            lines.append(f"  ΔS:          {p.ds:.2f} cal/(mol·K)")
            lines.append(f"  Score:       {p.score:.3f}")
            if p.warnings:
                lines.append(f"  Warnings:    {'; '.join(p.warnings)}")
            lines.append("")

        return "\n".join(lines)

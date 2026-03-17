"""Locked Nucleic Acid (LNA) incorporation module.

Implements best practices for LNA placement in capture probes:
  - LNA content: 20-50% of total probe length
  - No more than 3 consecutive LNA bases
  - Prefer LNA at mismatch-discriminating positions (G/C bases)
  - Avoid LNA in self-complementary regions
  - Place LNA in the central portion of the probe for stability
  - Each LNA base increases Tm by approximately 2-8°C

References:
    - Vester & Wengel (2004) Biochemistry 43:13233-13241
    - Kaur et al. (2006) Biochemistry 45:7347-7355
    - Koshkin et al. (1998) Tetrahedron 54:3607-3630
"""

from typing import Dict, List, Optional, Tuple

from probe_designer.core.utils import gc_content, self_complementarity_score

# Best practice constraints
MAX_CONSECUTIVE_LNA = 3
MIN_LNA_FRACTION = 0.20
MAX_LNA_FRACTION = 0.50
PREFERRED_LNA_BASES = {"G", "C"}  # LNA on G/C gives best Tm boost


class LNADesigner:
    """Design LNA modification patterns for capture probes.

    Determines optimal placement of LNA bases within a probe sequence
    following best practices for LNA incorporation.

    Args:
        min_lna_fraction: Minimum fraction of LNA bases (default 0.20).
        max_lna_fraction: Maximum fraction of LNA bases (default 0.50).
        max_consecutive: Maximum consecutive LNA bases (default 3).
        prefer_gc: Prefer placing LNA on G/C bases (default True).
    """

    def __init__(
        self,
        min_lna_fraction: float = MIN_LNA_FRACTION,
        max_lna_fraction: float = MAX_LNA_FRACTION,
        max_consecutive: int = MAX_CONSECUTIVE_LNA,
        prefer_gc: bool = True,
    ) -> None:
        self.min_lna_fraction = min_lna_fraction
        self.max_lna_fraction = max_lna_fraction
        self.max_consecutive = max_consecutive
        self.prefer_gc = prefer_gc

    def _score_position(self, sequence: str, position: int) -> float:
        """Score a position for LNA placement suitability.

        Higher scores indicate better positions for LNA incorporation.

        Args:
            sequence: Probe DNA sequence (uppercase).
            position: 0-based index in the sequence.

        Returns:
            Suitability score (higher is better).
        """
        n = len(sequence)
        base = sequence[position]
        score = 0.0

        # Prefer G/C bases (LNA on G/C has larger Tm increase)
        if self.prefer_gc and base in PREFERRED_LNA_BASES:
            score += 3.0

        # Prefer central positions over terminal ones
        # Central LNA bases contribute more to duplex stability
        distance_from_center = abs(position - n / 2.0) / (n / 2.0)
        score += 2.0 * (1.0 - distance_from_center)

        # Slight preference for evenly spaced positions
        # (every 2nd or 3rd base is a common LNA pattern)
        if position % 2 == 0:
            score += 0.5
        if position % 3 == 0:
            score += 0.3

        return score

    def _check_consecutive(
        self, positions: List[int], new_pos: int
    ) -> bool:
        """Check if adding new_pos would create too many consecutive LNAs.

        Args:
            positions: Currently selected LNA positions.
            new_pos: Candidate position to add.

        Returns:
            True if adding the position is safe.
        """
        test_positions = sorted(positions + [new_pos])
        consecutive = 1
        for i in range(1, len(test_positions)):
            if test_positions[i] == test_positions[i - 1] + 1:
                consecutive += 1
                if consecutive > self.max_consecutive:
                    return False
            else:
                consecutive = 1
        return True

    def design_lna_pattern(
        self,
        sequence: str,
        target_lna_count: Optional[int] = None,
    ) -> List[int]:
        """Determine optimal LNA positions for a probe sequence.

        Uses a greedy algorithm that scores each position and selects
        the best ones while respecting consecutive-LNA constraints.

        Args:
            sequence: Probe DNA sequence.
            target_lna_count: Desired number of LNA bases. If None,
                automatically calculated from sequence length and
                fraction constraints.

        Returns:
            Sorted list of 0-based LNA positions.
        """
        seq = sequence.upper()
        n = len(seq)

        # Determine target count
        if target_lna_count is None:
            # Aim for midpoint of allowed fraction range
            target_fraction = (self.min_lna_fraction + self.max_lna_fraction) / 2.0
            target_lna_count = max(1, round(n * target_fraction))

        # Clamp to valid range
        min_count = max(1, int(n * self.min_lna_fraction))
        max_count = min(n, int(n * self.max_lna_fraction))
        target_lna_count = max(min_count, min(max_count, target_lna_count))

        # Score all positions
        scored_positions = [
            (self._score_position(seq, i), i) for i in range(n)
        ]
        # Sort by score descending
        scored_positions.sort(key=lambda x: -x[0])

        # Greedily select positions
        selected: List[int] = []
        for _score, pos in scored_positions:
            if len(selected) >= target_lna_count:
                break
            if self._check_consecutive(selected, pos):
                selected.append(pos)

        return sorted(selected)

    def validate_lna_pattern(
        self, sequence: str, lna_positions: List[int]
    ) -> Dict[str, object]:
        """Validate an LNA modification pattern against best practices.

        Args:
            sequence: Probe DNA sequence.
            lna_positions: 0-based LNA position indices.

        Returns:
            Dictionary with validation results:
                - 'valid': bool
                - 'warnings': list of warning messages
                - 'lna_fraction': fraction of LNA bases
                - 'max_consecutive_lna': max consecutive LNA count
                - 'gc_lna_fraction': fraction of LNA bases on G/C
        """
        seq = sequence.upper()
        n = len(seq)
        warnings: List[str] = []

        if not lna_positions:
            return {
                "valid": True,
                "warnings": ["No LNA modifications specified"],
                "lna_fraction": 0.0,
                "max_consecutive_lna": 0,
                "gc_lna_fraction": 0.0,
            }

        # Check fraction
        lna_frac = len(lna_positions) / n
        if lna_frac < self.min_lna_fraction:
            warnings.append(
                f"LNA fraction ({lna_frac:.2f}) below minimum ({self.min_lna_fraction})"
            )
        if lna_frac > self.max_lna_fraction:
            warnings.append(
                f"LNA fraction ({lna_frac:.2f}) above maximum ({self.max_lna_fraction})"
            )

        # Check consecutive
        sorted_pos = sorted(lna_positions)
        max_consec = 1
        current_consec = 1
        for i in range(1, len(sorted_pos)):
            if sorted_pos[i] == sorted_pos[i - 1] + 1:
                current_consec += 1
                max_consec = max(max_consec, current_consec)
            else:
                current_consec = 1

        if max_consec > self.max_consecutive:
            warnings.append(
                f"Maximum consecutive LNA ({max_consec}) exceeds limit ({self.max_consecutive})"
            )

        # Check G/C preference
        gc_lna = sum(1 for p in lna_positions if seq[p] in PREFERRED_LNA_BASES)
        gc_lna_frac = gc_lna / len(lna_positions) if lna_positions else 0.0
        if gc_lna_frac < 0.4:
            warnings.append(
                f"Low GC-LNA fraction ({gc_lna_frac:.2f}); consider placing more LNA on G/C bases"
            )

        # Check positions within bounds
        for p in lna_positions:
            if p < 0 or p >= n:
                warnings.append(f"LNA position {p} is out of bounds (0-{n - 1})")

        valid = not any(
            "exceeds" in w or "out of bounds" in w for w in warnings
        )

        return {
            "valid": valid,
            "warnings": warnings,
            "lna_fraction": round(lna_frac, 3),
            "max_consecutive_lna": max_consec,
            "gc_lna_fraction": round(gc_lna_frac, 3),
        }

    def format_lna_sequence(
        self, sequence: str, lna_positions: List[int]
    ) -> str:
        """Format a sequence showing LNA positions with '+' prefix notation.

        Standard convention: LNA bases are prefixed with '+'.
        Example: 'A+CG+T' means C and T are LNA-modified.

        Args:
            sequence: DNA sequence.
            lna_positions: 0-based LNA positions.

        Returns:
            Formatted string with LNA notation.
        """
        lna_set = set(lna_positions)
        parts = []
        for i, base in enumerate(sequence.upper()):
            if i in lna_set:
                parts.append(f"+{base}")
            else:
                parts.append(base)
        return "".join(parts)

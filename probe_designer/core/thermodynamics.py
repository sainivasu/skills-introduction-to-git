"""Nearest-neighbor thermodynamic calculations for DNA and LNA-modified probes.

Implements the SantaLucia (1998) unified nearest-neighbor parameters for
DNA/DNA duplex thermodynamics, with corrections for LNA modifications.

References:
    - SantaLucia J Jr. (1998) Proc Natl Acad Sci USA 95:1460-1465
    - Owczarzy et al. (2004) Biochemistry 43:3537-3554
    - McTigue et al. (2004) Biochemistry 43:5388-5405 (LNA thermodynamics)
"""

import math
from typing import Dict, List, Optional, Tuple

# SantaLucia (1998) unified nearest-neighbor parameters for DNA/DNA duplexes.
# Format: {dinucleotide: (ΔH kcal/mol, ΔS cal/mol·K)}
DNA_NN_PARAMS: Dict[str, Tuple[float, float]] = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

# Initiation parameters for DNA duplexes
# (ΔH kcal/mol, ΔS cal/mol·K)
INIT_WITH_TERMINAL_GC = (0.1, -2.8)
INIT_WITH_TERMINAL_AT = (2.3, 4.1)

# Symmetry correction for self-complementary sequences
SYMMETRY_CORRECTION_S = -1.4  # cal/mol·K

# LNA modification thermodynamic bonus per LNA base (average values).
# Based on McTigue et al. (2004) and other studies.
# LNA bases increase ΔH (more negative = more stable) and ΔS (more negative).
LNA_DH_BONUS = -1.0   # kcal/mol per LNA base (additional stability)
LNA_DS_BONUS = -2.5   # cal/mol·K per LNA base

# Gas constant in cal/(mol·K)
R = 1.987


class ThermodynamicCalculator:
    """Calculate thermodynamic properties of DNA and LNA-modified probes.

    Uses nearest-neighbor parameters to compute melting temperature (Tm),
    Gibbs free energy (ΔG), enthalpy (ΔH), and entropy (ΔS) for
    oligonucleotide duplexes.

    Args:
        na_conc: Sodium ion concentration in molar (default 0.05 M).
        mg_conc: Magnesium ion concentration in molar (default 0.0 M).
        oligo_conc: Oligonucleotide concentration in molar (default 250e-9 M).
        target_conc: Target DNA concentration in molar (default 250e-9 M).
    """

    def __init__(
        self,
        na_conc: float = 0.05,
        mg_conc: float = 0.0,
        oligo_conc: float = 250e-9,
        target_conc: float = 250e-9,
    ) -> None:
        self.na_conc = na_conc
        self.mg_conc = mg_conc
        self.oligo_conc = oligo_conc
        self.target_conc = target_conc

    def _nn_parameters(
        self, sequence: str, lna_positions: Optional[List[int]] = None
    ) -> Tuple[float, float]:
        """Compute total ΔH and ΔS using nearest-neighbor model.

        Args:
            sequence: DNA sequence (uppercase).
            lna_positions: List of 0-based indices where LNA bases are placed.

        Returns:
            Tuple of (total_dH in kcal/mol, total_dS in cal/mol·K).
        """
        seq = sequence.upper()
        n = len(seq)

        if lna_positions is None:
            lna_positions = []
        lna_set = set(lna_positions)

        # Initiation parameters
        total_dh = 0.0
        total_ds = 0.0

        # 5' terminal initiation
        if seq[0] in ("G", "C"):
            total_dh += INIT_WITH_TERMINAL_GC[0]
            total_ds += INIT_WITH_TERMINAL_GC[1]
        else:
            total_dh += INIT_WITH_TERMINAL_AT[0]
            total_ds += INIT_WITH_TERMINAL_AT[1]

        # 3' terminal initiation
        if seq[-1] in ("G", "C"):
            total_dh += INIT_WITH_TERMINAL_GC[0]
            total_ds += INIT_WITH_TERMINAL_GC[1]
        else:
            total_dh += INIT_WITH_TERMINAL_AT[0]
            total_ds += INIT_WITH_TERMINAL_AT[1]

        # Nearest-neighbor stacking
        for i in range(n - 1):
            dinuc = seq[i : i + 2]
            if dinuc in DNA_NN_PARAMS:
                dh, ds = DNA_NN_PARAMS[dinuc]
                total_dh += dh
                total_ds += ds

        # LNA corrections
        num_lna = len(lna_set)
        if num_lna > 0:
            total_dh += num_lna * LNA_DH_BONUS
            total_ds += num_lna * LNA_DS_BONUS

        return total_dh, total_ds

    def _salt_correction(self, sequence_length: int) -> float:
        """Calculate salt correction factor for Tm using Owczarzy method.

        Returns a correction factor to add to 1/Tm.

        Args:
            sequence_length: Length of the oligonucleotide.

        Returns:
            Salt correction in 1/K units.
        """
        # Effective monovalent cation concentration
        # If Mg2+ is present, convert to equivalent Na+ contribution
        na_equiv = self.na_conc
        if self.mg_conc > 0:
            # Approximate: free Mg2+ contributes ~3.3x to ionic strength
            na_equiv += self.mg_conc * 3.3

        if na_equiv <= 0:
            na_equiv = 0.01  # minimum salt concentration

        # SantaLucia (1998) salt correction
        # 1/Tm(salt) = 1/Tm(1M) + (4.29 * f_GC - 3.95) * 1e-5 * ln[Na+]
        # + 9.40e-6 * (ln[Na+])^2
        # Simplified version: correction per base
        return (4.29 * 0.5 - 3.95) * 1e-5 * math.log(na_equiv)

    def calculate(
        self, sequence: str, lna_positions: Optional[List[int]] = None
    ) -> Dict[str, float]:
        """Calculate full thermodynamic profile for a probe-target duplex.

        Args:
            sequence: Probe DNA sequence (5' to 3').
            lna_positions: 0-based indices of LNA-modified positions.

        Returns:
            Dictionary with keys:
                - 'dH': Enthalpy in kcal/mol
                - 'dS': Entropy in cal/(mol·K)
                - 'dG_25': Gibbs free energy at 25°C in kcal/mol
                - 'dG_37': Gibbs free energy at 37°C in kcal/mol
                - 'Tm': Melting temperature in °C
                - 'num_lna': Number of LNA modifications
        """
        seq = sequence.upper()
        if lna_positions is None:
            lna_positions = []

        total_dh, total_ds = self._nn_parameters(seq, lna_positions)

        # Concentration correction for non-self-complementary
        # Ct = total strand concentration; for non-self-complementary:
        # if [probe] >> [target], Ct ~ [probe]; otherwise use
        # Ct = (C_probe - C_target/2) when probe is in excess
        ct = max(self.oligo_conc, self.target_conc)

        # Tm in Kelvin (nearest-neighbor formula)
        # Tm = ΔH / (ΔS + R * ln(Ct/4))
        ds_total = total_ds + R * math.log(ct / 4.0)

        if ds_total == 0:
            tm_kelvin = 0.0
        else:
            tm_kelvin = (total_dh * 1000.0) / ds_total

        # Apply salt correction
        if tm_kelvin > 0:
            salt_corr = self._salt_correction(len(seq))
            tm_kelvin_corr = 1.0 / (1.0 / tm_kelvin + salt_corr)
        else:
            tm_kelvin_corr = tm_kelvin

        tm_celsius = tm_kelvin_corr - 273.15

        # Free energies
        dg_25 = total_dh - (273.15 + 25.0) * (total_ds / 1000.0)
        dg_37 = total_dh - (273.15 + 37.0) * (total_ds / 1000.0)

        return {
            "dH": round(total_dh, 2),
            "dS": round(total_ds, 2),
            "dG_25": round(dg_25, 2),
            "dG_37": round(dg_37, 2),
            "Tm": round(tm_celsius, 1),
            "num_lna": len(lna_positions),
        }

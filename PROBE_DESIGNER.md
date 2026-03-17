# Hybrid Capture Probe Designer

A Python tool for designing ultra-short hybrid capture probes (<40 bp) with optional Locked Nucleic Acid (LNA) incorporation for enhanced binding affinity and specificity.

## Overview

This tool addresses the challenge of designing capture probes for ultra-short DNA fragments, where conventional probe design approaches may produce probes with insufficient melting temperatures (Tm) for stable hybridization. By incorporating LNA modifications following established best practices, the tool can boost probe Tm while maintaining high specificity.

### Key Features

- **Ultra-short probe design**: Probes from 15-40 bp targeting short DNA fragments
- **Nearest-neighbor thermodynamics**: SantaLucia (1998) unified parameters for accurate Tm, ΔG, ΔH, and ΔS calculations
- **LNA incorporation**: Automatic placement of LNA bases following best practices:
  - 20-50% LNA content for optimal binding enhancement
  - Maximum 3 consecutive LNA bases to avoid aggregation
  - Preferential placement on G/C bases for maximal Tm boost (~2-8°C per LNA)
  - Central positioning for improved duplex stability
- **Multi-criteria probe scoring**: Considers Tm optimality, GC content, thermodynamic stability, and self-complementarity
- **Both-strand design**: Generates probes from sense and antisense strands
- **Configurable parameters**: Salt concentration, probe concentration, GC bounds, Tm targets, and more

### Design Algorithm

1. **Tiling**: Generate candidate probes by sliding across the target sequence at configurable step sizes and multiple probe lengths
2. **Filtering**: Remove candidates with extreme GC content, homopolymer runs (>4), or high self-complementarity
3. **Thermodynamic evaluation**: Calculate Tm, ΔG, ΔH, ΔS using the nearest-neighbor model with salt corrections
4. **LNA incorporation** (optional): For probes with Tm below the target, add LNA modifications to boost stability
5. **Scoring and ranking**: Score probes on Tm proximity to target, GC content balance, thermodynamic stability, and structural quality

### References

- SantaLucia J Jr. (1998) "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." *PNAS* 95:1460-1465
- Owczarzy et al. (2004) "Effects of sodium ions on DNA duplex oligomers." *Biochemistry* 43:3537-3554
- McTigue et al. (2004) "Sequence-dependent thermodynamic parameters for locked nucleic acid." *Biochemistry* 43:5388-5405
- Vester & Wengel (2004) "LNA: a versatile tool for therapeutics and genomics." *Biochemistry* 43:13233-13241
- Kaur et al. (2006) "Thermodynamic, counterion, and hydration effects for LNA/DNA hybridization." *Biochemistry* 45:7347-7355

## Installation

```bash
# From the repository root:
pip install -e .
```

No external dependencies are required — the tool uses only the Python standard library (Python ≥ 3.8).

## Usage

### Command Line

```bash
# Basic usage — design probes for a 25 bp target
probe-designer ATCGATCGTAGCTAGCTTGACGTAC

# Or use python -m:
python -m probe_designer.cli ATCGATCGTAGCTAGCTTGACGTAC

# Return top 5 probes, disable LNA
probe-designer ATCGATCGTAGCTAGCTTGACGTAC -n 5 --no-lna

# Fixed probe length of 18 bp
probe-designer ATCGATCGTAGCTAGCTTGACGTAC -l 18

# Custom salt and Tm parameters
probe-designer ATCGATCGTAGCTAGCTTGACGTAC --na-conc 0.1 --optimal-tm 65

# From a FASTA file
probe-designer target.fasta -n 10
```

### CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `target` | (required) | Target DNA sequence or FASTA file path |
| `-l, --probe-length` | auto | Fixed probe length in bp |
| `-n, --top-n` | 10 | Number of top probes to return |
| `--min-length` | 15 | Minimum probe length (bp) |
| `--max-length` | 40 | Maximum probe length (bp) |
| `--min-gc` | 0.30 | Minimum GC content fraction |
| `--max-gc` | 0.70 | Maximum GC content fraction |
| `--min-tm` | 45.0 | Minimum melting temperature (°C) |
| `--max-tm` | 75.0 | Maximum melting temperature (°C) |
| `--optimal-tm` | 60.0 | Target melting temperature (°C) |
| `--na-conc` | 0.05 | Sodium concentration (M) |
| `--mg-conc` | 0.0 | Magnesium concentration (M) |
| `--no-lna` | False | Disable LNA incorporation |
| `--single-strand` | False | Design probes on sense strand only |
| `--tiling-step` | 1 | Tiling step size (bp) |

### Python API

```python
from probe_designer.core import ProbeDesigner, ThermodynamicCalculator, LNADesigner

# Design probes for a target
designer = ProbeDesigner(use_lna=True, optimal_tm=60.0)
probes = designer.design_probes("ATCGATCGTAGCTAGCTTGA", top_n=5)

for probe in probes:
    print(f"{probe.lna_sequence or probe.sequence}  Tm={probe.tm}°C  Score={probe.score}")

# Direct thermodynamic calculations
calc = ThermodynamicCalculator(na_conc=0.05)
result = calc.calculate("ATCGATCGTAGCTAGC", lna_positions=[2, 5, 8])
print(f"Tm = {result['Tm']}°C, ΔG(37°C) = {result['dG_37']} kcal/mol")

# LNA pattern design
lna = LNADesigner()
positions = lna.design_lna_pattern("ATCGATCGTAGCTAGC")
formatted = lna.format_lna_sequence("ATCGATCGTAGCTAGC", positions)
print(f"LNA pattern: {formatted}")
```

## Output Format

LNA-modified bases are shown with the `+` prefix notation (standard convention):

```
Sequence:    ATCGATCGTAGCTAGC
LNA format:  AT+C+GAT+C+GTA+G+CTAGC
```

In this example, positions marked with `+` are LNA-modified bases (`+C`, `+C`, `+G`, `+C`).

## Running Tests

```bash
python -m unittest discover -s tests -v
```

## Project Structure

```
probe_designer/
├── __init__.py              # Package metadata
├── cli.py                   # Command-line interface
└── core/
    ├── __init__.py          # Core module exports
    ├── designer.py          # Main probe design engine
    ├── lna.py               # LNA incorporation logic
    ├── thermodynamics.py    # Nearest-neighbor thermodynamic calculations
    └── utils.py             # Sequence utilities (GC, complement, validation)
tests/
├── test_cli.py              # CLI tests
├── test_designer.py         # Probe designer tests
├── test_lna.py              # LNA module tests
├── test_thermodynamics.py   # Thermodynamics tests
└── test_utils.py            # Utility function tests
setup.py                     # Package configuration
```

## LNA Best Practices Implemented

1. **Fraction control**: LNA content maintained at 20-50% of probe length
2. **Consecutive limit**: No more than 3 consecutive LNA bases (prevents aggregation and reduces synthesis issues)
3. **G/C preference**: LNA preferentially placed on G/C bases where Tm boost is greatest
4. **Central placement**: LNA bases positioned toward the probe center for maximal duplex stabilization
5. **Validation**: Built-in validation checks all patterns against best practice rules

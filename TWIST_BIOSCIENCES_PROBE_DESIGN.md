# Twist Biosciences Probe Design Technology

## Overview

Twist Biosciences uses a **silicon-based high-throughput DNA synthesis platform** to manufacture custom oligonucleotide probes for targeted hybridization capture. Their technology enables massively parallel synthesis of hundreds of thousands of unique probes on a miniaturized silicon chip, delivering high uniformity and accuracy at scale.

This document explains how Twist probes are designed to **hybridize specifically with complementary target strands** while **minimizing unwanted probe-to-probe (self-dimer) interactions**.

---

## How Probes Target Complementary Strands

### Hybridization Capture Principle

Twist probes are single-stranded DNA oligonucleotides (typically 80–120 nt) designed to be the **reverse complement** of a genomic region of interest. During a hybridization capture workflow:

1. **Genomic DNA is fragmented** and denatured into single strands.
2. **Biotinylated probes** are added and allowed to hybridize with their complementary target fragments under controlled temperature and salt conditions.
3. **Streptavidin-coated magnetic beads** pull down probe–target duplexes.
4. Non-target DNA is washed away, enriching for the regions of interest.

Specificity depends on **Watson–Crick base pairing** (A–T, G–C) governed by thermodynamic stability. Twist designs each probe so that its melting temperature (Tm) against the intended target is significantly higher than against any off-target sequence.

### Uniform Tiling and Coverage

Twist panels tile probes across every target region with overlapping coverage. This redundancy ensures that even if one probe performs sub-optimally, neighboring probes capture the same region.

---

## How Probes Avoid Binding to Each Other

A pool of thousands of probes risks **cross-hybridization** (one probe binding another) if their sequences share complementarity. Twist addresses this through multiple design strategies:

### 1. Computational Sequence Screening

During panel design, every candidate probe is screened against the entire probe pool for potential **self-complementarity** and **inter-probe complementarity**:

- **Self-dimer check** — A probe's sequence is aligned against its own reverse complement. Probes with long internal palindromes or hairpin-forming regions (typically ≥ 6 contiguous complementary bases) are flagged and redesigned.
- **Cross-dimer check** — Each probe is compared pairwise against all other probes in the pool. Pairs that share significant complementary stretches are identified. One or both probes are then shifted, shortened, or replaced with an alternative tiling position to break the complementarity.

### 2. Thermodynamic Balancing (Tm Optimization)

Probes are designed so that:

| Parameter | Target |
|---|---|
| **Probe–target Tm** | High (typically 60–75 °C), ensuring stable capture |
| **Probe–probe Tm** | Low (well below hybridization temperature), ensuring duplexes do not form |

By keeping the hybridization reaction temperature close to the probe–target Tm but far above any probe–probe Tm, unintended duplexes are thermodynamically disfavored and melt apart during the reaction.

### 3. GC Content and Complexity Filtering

- Probes with **extreme GC content** (< 30% or > 70%) are avoided because they either bind too weakly (AT-rich) or form excessively stable secondary structures (GC-rich).
- **Low-complexity sequences** (e.g., homopolymer runs, dinucleotide repeats) are excluded because they are prone to non-specific binding and cross-hybridization.

### 4. Blocking Agents

In the wet-lab protocol, **Cot-1 DNA** (enriched for repetitive human sequences) and **blocking oligonucleotides** complementary to adapter sequences are added to the hybridization reaction. These blockers competitively bind repetitive and adapter-derived sequences, preventing probes from interacting with each other through shared adapter tails or repetitive content.

### 5. Silicon-Platform Synthesis Precision

Twist's proprietary silicon synthesis platform writes each probe independently in a microwell on the chip surface, achieving:

- **High sequence fidelity** (low error rate ~1 in 2000 bases), so probes match their designed sequence and behave as computationally predicted.
- **Uniform representation** across the pool, preventing any single probe from being over-represented and driving off-target interactions by mass action.

---

## Summary

| Challenge | Twist Solution |
|---|---|
| Target specificity | Reverse-complement design with Tm-optimized hybridization |
| Self-dimers / hairpins | Computational screening for internal complementarity |
| Cross-hybridization between probes | Pairwise complementarity checks and sequence redesign |
| Repetitive sequence interference | Cot-1 DNA and adapter blockers in the reaction |
| Probe quality and uniformity | Silicon-based high-fidelity parallel synthesis |

By combining **computational probe design** (sequence screening, thermodynamic modeling, complexity filtering) with **wet-lab blocking strategies** and **high-fidelity silicon synthesis**, Twist Biosciences ensures that probes bind their intended genomic targets with high specificity while avoiding unproductive probe-to-probe interactions.

---

## References

- Twist Biosciences. *Target Enrichment with Twist Custom Panels.* [twist.com](https://www.twistbioscience.com)
- SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS*, 95(4), 1460–1465.
- Gnirke, A. et al. (2009). Solution hybrid selection with ultra-long oligonucleotides for massively parallel targeted sequencing. *Nature Biotechnology*, 27(2), 182–189.

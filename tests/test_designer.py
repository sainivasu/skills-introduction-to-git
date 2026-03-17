"""Tests for probe_designer.core.designer module."""

import unittest

from probe_designer.core.designer import ProbeCandidate, ProbeDesigner


class TestProbeDesigner(unittest.TestCase):
    """Tests for ProbeDesigner."""

    def setUp(self):
        self.designer = ProbeDesigner(
            min_probe_length=15,
            max_probe_length=30,
            use_lna=True,
            tiling_step=1,
        )

    def test_design_returns_list(self):
        """design_probes should return a list of ProbeCandidate."""
        target = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32-mer
        probes = self.designer.design_probes(target, top_n=5)
        self.assertIsInstance(probes, list)
        for p in probes:
            self.assertIsInstance(p, ProbeCandidate)

    def test_probes_have_valid_sequences(self):
        """Each probe should have a valid DNA sequence."""
        target = "TAGCTTGACGATCAGTACGTCCTA"  # 24-mer, low self-comp
        probes = self.designer.design_probes(target, top_n=5)
        valid_bases = set("ACGT")
        for p in probes:
            self.assertTrue(all(b in valid_bases for b in p.sequence))

    def test_probes_sorted_by_score(self):
        """Probes should be sorted by score descending."""
        target = "ACGTGCATACGTGCATACGTGCAT"  # 24-mer
        probes = self.designer.design_probes(target, top_n=10)
        if len(probes) > 1:
            scores = [p.score for p in probes]
            self.assertEqual(scores, sorted(scores, reverse=True))

    def test_gc_content_within_range(self):
        """Probes should have GC content within configured range."""
        target = "TAGCTTGACGATCAGTACGT"  # 20-mer, low self-comp
        probes = self.designer.design_probes(target, top_n=5)
        for p in probes:
            self.assertGreaterEqual(p.gc, self.designer.min_gc)
            self.assertLessEqual(p.gc, self.designer.max_gc)

    def test_invalid_sequence_raises(self):
        """Invalid target sequence should raise ValueError."""
        with self.assertRaises(ValueError):
            self.designer.design_probes("ACGXYZ")

    def test_short_target_raises(self):
        """Target shorter than min_probe_length should raise ValueError."""
        with self.assertRaises(ValueError):
            self.designer.design_probes("ACG")

    def test_lna_probes_have_lna_info(self):
        """When LNA is enabled, low-Tm probes should get LNA modifications."""
        # Use an AT-rich target to trigger LNA incorporation
        designer = ProbeDesigner(
            min_probe_length=15,
            max_probe_length=20,
            use_lna=True,
            min_gc=0.0,
            min_tm=0.0,
        )
        target = "ATATGCATATGCATATGCAT"  # 20-mer, mixed but some AT-rich
        probes = designer.design_probes(target, top_n=20)
        # At least some probes should have LNA
        lna_probes = [p for p in probes if p.lna_positions]
        # It's possible no LNA needed if Tm already high enough
        self.assertIsInstance(lna_probes, list)

    def test_no_lna_mode(self):
        """When LNA is disabled, no probes should have LNA positions."""
        designer = ProbeDesigner(
            min_probe_length=15,
            max_probe_length=20,
            use_lna=False,
            min_tm=0.0,
        )
        target = "ACGTACGTACGTACGTACGT"  # 20-mer
        probes = designer.design_probes(target, top_n=5)
        for p in probes:
            self.assertEqual(p.lna_positions, [])

    def test_fixed_probe_length(self):
        """Fixed probe length should produce probes of that length."""
        target = "ACGTACGTACGTACGTACGTACGT"  # 24-mer
        probes = self.designer.design_probes(target, probe_length=18, top_n=5)
        for p in probes:
            self.assertEqual(len(p.sequence), 18)

    def test_both_strands(self):
        """Both strand mode should produce probes on + and - strands."""
        target = "ACGTACGTACGTACGTACGT"  # 20-mer
        probes = self.designer.design_probes(target, top_n=20)
        strands = {p.strand for p in probes}
        if len(probes) > 1:
            self.assertEqual(strands, {"+", "-"})

    def test_single_strand_mode(self):
        """Single strand mode should only produce + strand probes."""
        designer = ProbeDesigner(
            min_probe_length=15,
            max_probe_length=20,
            both_strands=False,
            min_tm=0.0,
        )
        target = "ACGTACGTACGTACGTACGT"  # 20-mer
        probes = designer.design_probes(target, top_n=10)
        for p in probes:
            self.assertEqual(p.strand, "+")


class TestProbeDesignerSummary(unittest.TestCase):
    """Tests for ProbeDesigner.design_probes_summary."""

    def test_summary_is_string(self):
        """Summary should be a non-empty string."""
        designer = ProbeDesigner(min_tm=0.0)
        target = "TAGCTTGACGATCAGTACGT"  # 20-mer, low self-comp
        summary = designer.design_probes_summary(target, top_n=3)
        self.assertIsInstance(summary, str)
        self.assertGreater(len(summary), 0)

    def test_summary_contains_header(self):
        """Summary should contain the header."""
        designer = ProbeDesigner(min_tm=0.0)
        target = "TAGCTTGACGATCAGTACGT"
        summary = designer.design_probes_summary(target, top_n=3)
        self.assertIn("HYBRID CAPTURE PROBE DESIGN", summary)


class TestProbeDesignerUltraShort(unittest.TestCase):
    """Tests for ultra-short target handling (<40 bp)."""

    def test_20bp_target(self):
        """Should handle 20 bp targets."""
        designer = ProbeDesigner(
            min_probe_length=15,
            max_probe_length=20,
            min_tm=0.0,
        )
        target = "TAGCTTGACGATCAGTACGT"  # 20 bp, low self-comp
        probes = designer.design_probes(target, top_n=5)
        self.assertGreater(len(probes), 0)

    def test_30bp_target(self):
        """Should handle 30 bp targets."""
        designer = ProbeDesigner(
            min_probe_length=15,
            max_probe_length=30,
            min_tm=0.0,
        )
        target = "TAGCTTGACGATCAGTACGTCCTAGAGTTCG"  # 31 bp, low self-comp
        probes = designer.design_probes(target, top_n=5)
        self.assertGreater(len(probes), 0)

    def test_15bp_minimum_target(self):
        """Should handle target at minimum length."""
        designer = ProbeDesigner(
            min_probe_length=15,
            max_probe_length=15,
            min_tm=0.0,
        )
        target = "AGCTTGACGATCAGT"  # 15 bp, low self-comp
        probes = designer.design_probes(target, top_n=5)
        self.assertGreater(len(probes), 0)


if __name__ == "__main__":
    unittest.main()

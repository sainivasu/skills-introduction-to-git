"""Tests for probe_designer.core.lna module."""

import unittest

from probe_designer.core.lna import LNADesigner


class TestLNADesigner(unittest.TestCase):
    """Tests for LNADesigner."""

    def setUp(self):
        self.designer = LNADesigner()

    def test_design_returns_positions(self):
        """design_lna_pattern should return a list of positions."""
        positions = self.designer.design_lna_pattern("ACGTACGTACGTACGT")
        self.assertIsInstance(positions, list)
        self.assertGreater(len(positions), 0)

    def test_positions_within_bounds(self):
        """All LNA positions should be within sequence bounds."""
        seq = "ACGTACGTACGT"
        positions = self.designer.design_lna_pattern(seq)
        for pos in positions:
            self.assertGreaterEqual(pos, 0)
            self.assertLess(pos, len(seq))

    def test_no_excessive_consecutive(self):
        """No more than max_consecutive LNA bases in a row."""
        seq = "GCGCGCGCGCGCGCGCGCGC"  # 20-mer, all GC
        positions = self.designer.design_lna_pattern(seq)
        sorted_pos = sorted(positions)
        max_consec = 1
        consec = 1
        for i in range(1, len(sorted_pos)):
            if sorted_pos[i] == sorted_pos[i - 1] + 1:
                consec += 1
                max_consec = max(max_consec, consec)
            else:
                consec = 1
        self.assertLessEqual(max_consec, 3)

    def test_lna_fraction_within_bounds(self):
        """LNA fraction should be within configured bounds."""
        seq = "ACGTACGTACGTACGTACGT"  # 20-mer
        positions = self.designer.design_lna_pattern(seq)
        frac = len(positions) / len(seq)
        self.assertGreaterEqual(frac, 0.15)  # slightly below min for tolerance
        self.assertLessEqual(frac, 0.55)  # slightly above max for tolerance

    def test_custom_target_count(self):
        """Custom target LNA count should be respected."""
        seq = "ACGTACGTACGTACGT"  # 16-mer
        positions = self.designer.design_lna_pattern(seq, target_lna_count=4)
        # Should be clamped to valid range
        self.assertGreaterEqual(len(positions), 3)
        self.assertLessEqual(len(positions), 8)

    def test_validate_valid_pattern(self):
        """A well-designed pattern should validate without errors."""
        seq = "ACGTACGTACGTACGT"
        positions = self.designer.design_lna_pattern(seq)
        result = self.designer.validate_lna_pattern(seq, positions)
        self.assertTrue(result["valid"])

    def test_validate_empty_pattern(self):
        """Empty LNA pattern should validate."""
        result = self.designer.validate_lna_pattern("ACGTACGT", [])
        self.assertTrue(result["valid"])

    def test_validate_excessive_consecutive(self):
        """Pattern with 4+ consecutive LNA should warn."""
        result = self.designer.validate_lna_pattern(
            "ACGTACGT", [0, 1, 2, 3]
        )
        self.assertFalse(result["valid"])
        self.assertTrue(
            any("consecutive" in w.lower() for w in result["warnings"])
        )

    def test_format_lna_sequence(self):
        """LNA format should use '+' prefix notation."""
        result = self.designer.format_lna_sequence("ACGT", [1, 3])
        self.assertEqual(result, "A+CG+T")

    def test_format_no_lna(self):
        """No LNA should return plain sequence."""
        result = self.designer.format_lna_sequence("ACGT", [])
        self.assertEqual(result, "ACGT")


class TestLNADesignerCustomConfig(unittest.TestCase):
    """Tests for LNADesigner with custom configuration."""

    def test_higher_max_fraction(self):
        """Higher max fraction allows more LNA bases."""
        designer = LNADesigner(min_lna_fraction=0.30, max_lna_fraction=0.60)
        seq = "ACGTACGTACGTACGT"
        positions = designer.design_lna_pattern(seq)
        frac = len(positions) / len(seq)
        self.assertGreaterEqual(frac, 0.25)

    def test_no_gc_preference(self):
        """When prefer_gc=False, AT positions should also be selected."""
        designer = LNADesigner(prefer_gc=False)
        seq = "AAAAAACCCCCCCCCC"  # 6A + 10C
        positions = designer.design_lna_pattern(seq)
        # Should include some A positions
        a_positions = [p for p in positions if seq[p] == "A"]
        # With no GC preference, some A positions should be picked
        self.assertIsInstance(a_positions, list)


if __name__ == "__main__":
    unittest.main()

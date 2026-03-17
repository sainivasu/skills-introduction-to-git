"""Tests for probe_designer.core.thermodynamics module."""

import unittest

from probe_designer.core.thermodynamics import ThermodynamicCalculator


class TestThermodynamicCalculator(unittest.TestCase):
    """Tests for ThermodynamicCalculator."""

    def setUp(self):
        self.calc = ThermodynamicCalculator()

    def test_basic_calculation(self):
        """Test that calculation returns expected keys."""
        result = self.calc.calculate("ACGTACGTAC")
        self.assertIn("dH", result)
        self.assertIn("dS", result)
        self.assertIn("dG_25", result)
        self.assertIn("dG_37", result)
        self.assertIn("Tm", result)
        self.assertIn("num_lna", result)

    def test_tm_in_reasonable_range(self):
        """Tm for a 20-mer should be in a reasonable range."""
        result = self.calc.calculate("ACGTACGTACGTACGTACGT")
        self.assertTrue(20 < result["Tm"] < 80)

    def test_gc_rich_higher_tm(self):
        """GC-rich sequences should have higher Tm."""
        gc_rich = self.calc.calculate("GCGCGCGCGCGCGCGCGCGC")
        at_rich = self.calc.calculate("ATATATATATATATATATATAT")
        self.assertGreater(gc_rich["Tm"], at_rich["Tm"])

    def test_dg_negative(self):
        """Free energy should be negative for stable duplexes."""
        result = self.calc.calculate("GCGCGCGCGC")
        self.assertLess(result["dG_37"], 0)

    def test_lna_increases_tm(self):
        """LNA modifications should increase Tm."""
        seq = "ACGTACGTACGTACGT"
        result_dna = self.calc.calculate(seq)
        result_lna = self.calc.calculate(seq, lna_positions=[2, 5, 8, 11])
        self.assertGreater(result_lna["Tm"], result_dna["Tm"])

    def test_lna_count_reported(self):
        """num_lna should match the number of LNA positions provided."""
        result = self.calc.calculate("ACGTACGT", lna_positions=[1, 3, 5])
        self.assertEqual(result["num_lna"], 3)

    def test_no_lna_count_zero(self):
        """num_lna should be 0 when no LNA positions provided."""
        result = self.calc.calculate("ACGTACGT")
        self.assertEqual(result["num_lna"], 0)

    def test_longer_sequence_lower_dg(self):
        """Longer sequences should have more negative ΔG (more stable)."""
        short = self.calc.calculate("GCGCGC")
        long = self.calc.calculate("GCGCGCGCGCGCGCGC")
        self.assertLess(long["dG_37"], short["dG_37"])

    def test_custom_salt_concentration(self):
        """Different salt concentrations should affect Tm."""
        calc_low = ThermodynamicCalculator(na_conc=0.01)
        calc_high = ThermodynamicCalculator(na_conc=0.5)
        seq = "ACGTACGTACGTACGT"
        tm_low = calc_low.calculate(seq)["Tm"]
        tm_high = calc_high.calculate(seq)["Tm"]
        # Higher salt generally increases Tm
        self.assertNotEqual(tm_low, tm_high)


if __name__ == "__main__":
    unittest.main()

"""Tests for probe_designer.core.utils module."""

import unittest

from probe_designer.core.utils import (
    gc_content,
    has_homopolymer_run,
    melt_temp_basic,
    reverse_complement,
    self_complementarity_score,
    validate_sequence,
)


class TestValidateSequence(unittest.TestCase):
    """Tests for validate_sequence."""

    def test_valid_uppercase(self):
        self.assertTrue(validate_sequence("ACGT"))

    def test_valid_lowercase(self):
        self.assertTrue(validate_sequence("acgt"))

    def test_valid_mixed_case(self):
        self.assertTrue(validate_sequence("AcGt"))

    def test_invalid_characters(self):
        self.assertFalse(validate_sequence("ACGN"))
        self.assertFalse(validate_sequence("ACGU"))

    def test_empty_string(self):
        self.assertFalse(validate_sequence(""))


class TestReverseComplement(unittest.TestCase):
    """Tests for reverse_complement."""

    def test_basic(self):
        self.assertEqual(reverse_complement("ACGT"), "ACGT")

    def test_poly_a(self):
        self.assertEqual(reverse_complement("AAAA"), "TTTT")

    def test_asymmetric(self):
        self.assertEqual(reverse_complement("AACG"), "CGTT")

    def test_case_insensitive(self):
        result = reverse_complement("acgt")
        self.assertEqual(result.upper(), "ACGT")

    def test_invalid_raises(self):
        with self.assertRaises(ValueError):
            reverse_complement("ACGN")


class TestGCContent(unittest.TestCase):
    """Tests for gc_content."""

    def test_all_gc(self):
        self.assertAlmostEqual(gc_content("GCGC"), 1.0)

    def test_all_at(self):
        self.assertAlmostEqual(gc_content("ATAT"), 0.0)

    def test_50_percent(self):
        self.assertAlmostEqual(gc_content("ACGT"), 0.5)

    def test_empty_raises(self):
        with self.assertRaises(ValueError):
            gc_content("")


class TestMeltTempBasic(unittest.TestCase):
    """Tests for melt_temp_basic (Wallace rule)."""

    def test_short_oligo(self):
        # 10-mer: ACGTACGTAC -> 5 GC, 5 AT
        # Tm = 2*5 + 4*5 = 30
        tm = melt_temp_basic("ACGTACGTAC")
        self.assertAlmostEqual(tm, 30.0)

    def test_longer_oligo(self):
        seq = "ACGTACGTACGTACGT"  # 16-mer, 8 GC, 8 AT
        tm = melt_temp_basic(seq)
        self.assertIsInstance(tm, float)
        # Should use the longer formula
        self.assertTrue(20 < tm < 80)


class TestSelfComplementarityScore(unittest.TestCase):
    """Tests for self_complementarity_score."""

    def test_palindrome(self):
        # ACGT is its own reverse complement
        score = self_complementarity_score("ACGT")
        self.assertGreater(score, 0)

    def test_poly_a(self):
        # AAAA vs TTTT - no self-complementarity at same offset
        score = self_complementarity_score("AAAA")
        self.assertIsInstance(score, int)

    def test_no_self_comp(self):
        # All same base
        score = self_complementarity_score("GGGG")
        # RC is CCCC, so no matches
        self.assertEqual(score, 0)


class TestHasHomopolymerRun(unittest.TestCase):
    """Tests for has_homopolymer_run."""

    def test_no_run(self):
        self.assertFalse(has_homopolymer_run("ACGTACGT", max_run=4))

    def test_has_run(self):
        self.assertTrue(has_homopolymer_run("AAAAACGT", max_run=4))

    def test_exact_boundary(self):
        self.assertFalse(has_homopolymer_run("AAAACGT", max_run=4))
        self.assertTrue(has_homopolymer_run("AAAAACGT", max_run=4))


if __name__ == "__main__":
    unittest.main()

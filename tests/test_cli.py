"""Tests for probe_designer.cli module."""

import unittest

from probe_designer.cli import main, parse_args


class TestParseArgs(unittest.TestCase):
    """Tests for CLI argument parsing."""

    def test_basic_target(self):
        args = parse_args(["ACGTACGTACGTACGT"])
        self.assertEqual(args.target, "ACGTACGTACGTACGT")

    def test_probe_length(self):
        args = parse_args(["ACGT", "-l", "20"])
        self.assertEqual(args.probe_length, 20)

    def test_top_n(self):
        args = parse_args(["ACGT", "-n", "5"])
        self.assertEqual(args.top_n, 5)

    def test_no_lna_flag(self):
        args = parse_args(["ACGT", "--no-lna"])
        self.assertTrue(args.no_lna)

    def test_defaults(self):
        args = parse_args(["ACGT"])
        self.assertIsNone(args.probe_length)
        self.assertEqual(args.top_n, 10)
        self.assertFalse(args.no_lna)
        self.assertFalse(args.single_strand)


class TestMainFunction(unittest.TestCase):
    """Tests for CLI main function."""

    def test_valid_sequence(self):
        """Main should return 0 for valid input."""
        result = main(["TAGCTTGACGATCAGTACGT", "-n", "3"])
        self.assertEqual(result, 0)

    def test_invalid_sequence(self):
        """Main should return 1 for invalid input."""
        result = main(["XYZXYZXYZ"])
        self.assertEqual(result, 1)

    def test_short_sequence(self):
        """Main should return 1 for too-short input."""
        result = main(["ACG"])
        self.assertEqual(result, 1)

    def test_no_lna_flag(self):
        """Main should work with --no-lna flag."""
        result = main(["TAGCTTGACGATCAGTACGT", "--no-lna", "-n", "3"])
        self.assertEqual(result, 0)

    def test_missing_file(self):
        """Main should return 1 for non-existent FASTA file."""
        result = main(["/nonexistent/file.fasta"])
        self.assertEqual(result, 1)


if __name__ == "__main__":
    unittest.main()

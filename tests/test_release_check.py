"""Publication regression checks use synthetic strings, never source records."""

from pathlib import Path
import tempfile
import unittest

from scripts.check_release import content_issues, inspect_tree


class ReleaseCheckTests(unittest.TestCase):
    def test_small_public_reference_tree_passes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "README.md").write_text("Public reference analysis.\n", encoding="utf-8")
            data = root / "data" / "public"
            data.mkdir(parents=True)
            (data / "evidence.csv").write_text("gene,source\nBUB1B,public literature\n", encoding="utf-8")
            count, issues = inspect_tree(root)
            self.assertEqual(count, 2)
            self.assertEqual(issues, [])

    def test_rejects_ignored_raw_files_and_archives(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / ".gitignore").write_text("*.gz\n", encoding="utf-8")
            (root / "source.vcf.gz").write_bytes(b"synthetic")
            _, issues = inspect_tree(root)
            self.assertTrue(any("source.vcf.gz" in issue for issue in issues))

    def test_renaming_raw_content_does_not_hide_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / "data" / "public"
            target.mkdir(parents=True)
            (target / "renamed.csv").write_text("##fileformat=" + "VCFv4.2\n", encoding="utf-8")
            _, issues = inspect_tree(root)
            self.assertTrue(any("VCF content" in issue for issue in issues))

    def test_secrets_report_category_without_echoing_value(self):
        synthetic_token = "gh" + "p_" + "a" * 36
        self.assertEqual(content_issues(synthetic_token), ["GitHub token"])
        self.assertNotIn(synthetic_token, " ".join(content_issues(synthetic_token)))

    def test_unknown_binary_and_oversize_files_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "unexpected.txt").write_text("unknown", encoding="utf-8")
            (root / "README.md").write_bytes(b"\xff\xfe")
            (root / "LICENSE").write_text("x" * 2_000_001, encoding="utf-8")
            _, issues = inspect_tree(root)
            self.assertEqual(len(issues), 3)

    def test_symlink_to_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / "README.md"
            target.write_text("Public data", encoding="utf-8")
            try:
                (root / "LICENSE").symlink_to(target)
            except (OSError, NotImplementedError):
                self.skipTest("Creating symlinks unavailable on this platform")
            _, issues = inspect_tree(root)
            self.assertTrue(any("symbolic links" in issue for issue in issues))


if __name__ == "__main__":
    unittest.main()

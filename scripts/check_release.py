"""Fail closed on unexpected files before publishing this small research repo.

Checks the working tree, including ignored data files, excluding Git metadata
and named development environments/caches. Do not force-add those excluded
directories to Git. This is not a history scanner.
This is a guardrail for accidental inclusions, not a substitute for reviewing
provenance and consent. It cannot establish that arbitrary text is deidentified.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys


ROOT_FILES = {
    ".gitignore", ".gitattributes", "README.md", "LICENSE", "LICENSE.md",
    "NOTICE", "CITATION.cff", "pyproject.toml", "requirements.txt",
    "requirements-dev.txt", "CHANGELOG.md", "CONTRIBUTING.md", "SECURITY.md",
}
TEXT_SUFFIXES = {".md", ".csv", ".tsv", ".json", ".svg"}
SKIP_DIRECTORIES = {
    ".git", ".venv", "venv", "__pycache__", ".pytest_cache", ".mypy_cache",
    ".ruff_cache", ".idea", ".vscode",
}
FORBIDDEN_DIRECTORIES = {
    "raw", "private", "clinical", "patient_data", "refs", "references",
    "downloads",
}
FORBIDDEN_SUFFIXES = {
    ".fastq", ".fq", ".fasta", ".fa", ".fna", ".faa", ".sam", ".bam",
    ".bai", ".cram", ".crai", ".vcf", ".bcf", ".gvcf", ".bed", ".bim",
    ".fam", ".pgen", ".pvar", ".psam", ".gz", ".bgz", ".zip", ".7z",
    ".tar", ".pdf", ".doc", ".docx", ".xlsx", ".xls", ".pem", ".key",
    ".p12", ".pfx",
}
SECRET_PATTERNS = (
    ("private key", re.compile(r"-----BEGIN (?:[A-Z ]+ )?PRIVATE KEY-----")),
    ("GitHub token", re.compile(r"\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,})\b")),
    ("AWS access key", re.compile(r"\b(?:AKIA|ASIA)[A-Z0-9]{16}\b")),
    ("Slack token", re.compile(r"\bxox[baprs]-[A-Za-z0-9-]{20,}\b")),
    ("API secret", re.compile(r"\bsk-(?:proj-)?[A-Za-z0-9_-]{30,}\b")),
    ("credential URL", re.compile(r"https?://[^\s/:@]+:[^\s/@]{8,}@")),
    ("literal credential", re.compile(
        r'''(?i)["']?(?:password|api_key|access_token|client_secret)["']?\s*[:=]\s*["'][^"'\s]{12,}["']'''
    )),
)
MAX_FILE_BYTES = 2_000_000


def allowed_path(relative: Path) -> bool:
    """Only source, documentation, and explicitly public small text data."""
    parts = relative.parts
    if len(parts) == 1:
        return relative.name in ROOT_FILES
    if parts[0] in {"src", "scripts", "tests"}:
        return relative.suffix == ".py"
    if parts[0] == "docs":
        return relative.suffix in {".md", ".svg"}
    if relative.as_posix() == "data/scenarios.json":
        return True
    if parts[:2] == ("data", "public") or parts[0] == "results":
        return relative.suffix in TEXT_SUFFIXES
    if parts[:2] == (".github", "workflows") and len(parts) == 3:
        return relative.suffix in {".yml", ".yaml"}
    return False


def content_issues(content: str) -> list[str]:
    issues = [label for label, pattern in SECRET_PATTERNS if pattern.search(content)]
    # Raw format headers remain detectable even if the file was renamed .csv.
    if re.search(r"(?m)^##fileformat=VCF", content):
        issues.append("VCF content")
    if re.search(r"(?m)^@HD\tVN:", content) or re.search(r"(?m)^@SQ\tSN:", content):
        issues.append("SAM content")
    if re.search(r"(?m)^@[^\n]+\n[ACGTNacgtn]{20,}\n\+[^\n]*\n[^\n]+", content):
        issues.append("FASTQ content")
    if re.search(r"(?m)^>[^\n]+\n[ACGTNacgtn]{40,}(?:\n|$)", content):
        issues.append("FASTA content")
    if "\x00" in content:
        issues.append("binary content")
    return issues


def inspect_tree(root: Path) -> tuple[int, list[str]]:
    issues: list[str] = []
    checked = 0

    def walk(directory: Path) -> None:
        nonlocal checked
        try:
            entries = sorted(directory.iterdir())
        except OSError as error:
            issues.append(f"{directory.relative_to(root).as_posix()}: cannot read directory ({error.__class__.__name__})")
            return
        for path in entries:
            relative = path.relative_to(root)
            label = relative.as_posix()
            # Never follow a link, including one with an otherwise excluded name.
            if path.is_symlink() or path.is_junction():
                issues.append(f"{label}: symbolic links and junctions are not permitted")
                continue
            if path.is_dir():
                if path.name in SKIP_DIRECTORIES or path.name.endswith(".egg-info"):
                    continue
                if path.name.lower() in FORBIDDEN_DIRECTORIES:
                    issues.append(f"{label}: source/private data directory is not permitted")
                    continue
                walk(path)
                continue
            checked += 1
            if not path.is_file():
                issues.append(f"{label}: unexpected filesystem object")
                continue
            if any(suffix.lower() in FORBIDDEN_SUFFIXES for suffix in path.suffixes):
                issues.append(f"{label}: source data, archive, clinical document, or key format")
                continue
            if not allowed_path(relative):
                issues.append(f"{label}: not on the publication allowlist")
                continue
            try:
                if path.stat().st_size > MAX_FILE_BYTES:
                    issues.append(f"{label}: exceeds {MAX_FILE_BYTES:,}-byte size limit")
                    continue
                content = path.read_text(encoding="utf-8-sig")
            except (OSError, UnicodeError) as error:
                issues.append(f"{label}: must be readable UTF-8 text ({error.__class__.__name__})")
                continue
            for issue in content_issues(content):
                # Report categories only; never echo a matched secret or source record.
                issues.append(f"{label}: {issue}")

    walk(root)
    return checked, issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        parser.error("--root must be an existing directory")
    checked, issues = inspect_tree(root)
    if issues:
        print(f"Release check FAILED ({len(issues)} issue(s), {checked} file(s) inspected).", file=sys.stderr)
        for issue in issues:
            print(f"- {issue}", file=sys.stderr)
        return 1
    print(f"Release check passed: {checked} allowed text files inspected.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

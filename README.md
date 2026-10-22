# BUB1B rescue map

Which molecular rescue approaches deserve testing for a BUB1B stop/missense model of mosaic variegated aneuploidy?

This repository turns that question into a reproducible sequence analysis and a research decision plan. It expands an earlier ataluren nomination to include native-protein correction, RNA and gene replacement, suppressor tRNAs, and exon skipping. It does not predict that any intervention will help a particular person.

**Track 1 is settled: compound-heterozygous BUB1B is the accepted diagnosis.** The rescue analysis carries forward **NM_001211.6:c.2210T>G (L737*) and c.3006T>G (N1002K), in trans**, from the completed diagnostic work. The participant reconfirmed that work on 22 October 2026 and supplied the Track 1 leaderboard context. Diagnosis is an input to this project; the open questions concern functional rescue. This repository contains public reference sequences and the accepted allele pair, without source genomic data or clinical records.

The [official Track 1 leaderboard](https://sagebio-rare-disease-real-kid-mva-hackathon-2026.hf.space/) was checked on 22 October 2026 and contains the supplied full-score BUB1B entries. Leaderboard performance concerns variant identification; therapeutic rescue is a separate research question.

## Start here

- [Research priorities and the decisive comparisons](docs/RESEARCH_PLAN.md)
- [Computed sequence results](results/rescue_map.md) and [machine-readable results](results/rescue_map.json)
- [Primary evidence and its limits](docs/EVIDENCE.md)
- [Data scope and repository status](docs/DATA_POLICY.md)

## What the analysis establishes

1. Standard ABE/CBE editing cannot recover either native residue in these codons. Exact reversions require G>T; N1002K also permits G>C to recover Asn with an alternative synonymous codon. Prime editing and newer transversion editors need separate feasibility work.
2. Conventional ADAR editing of the L737* UGA makes UGG, encoding **tryptophan**, not native leucine. Stop removal and restoration of normal protein are distinct outcomes.
3. Precisely skipping exon 17 removes 141 coding bases and produces a 1,003-residue protein. Explicit junction translation verifies deletion of residues **Glu715 through Leu761**. Function is unknown. Nonselective skipping also changes the N1002K copy and could impair its residual function.
4. Six exon junctions lie downstream of the stop. The positional nonsense-mediated decay criterion is met; the amount of surviving RNA is unmeasured.

These are sequence consequences, not measured editing efficiency, functional rescue or clinical benefit. The highest-value next evidence is whether defined protein products restore accurate chromosome segregation in the relevant genetic background.

## Reproduce locally

Python 3.12 or newer; no third-party packages needed. The checked-in reference snapshot permits an offline run.

```sh
python src/rescue_map.py
python -m unittest discover -s tests -v
python scripts/check_release.py
```

To rebuild the public reference snapshot from MANE 1.5, RefSeq NM_001211.6 and UniProt O60566:

```sh
python scripts/build_public_reference.py
```

The builder also accepts `--gtf /path/to/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`. It checks the full RefSeq translation against UniProt and records source hashes. Upstream annotation metadata can change; the committed snapshot and input hashes define the reproduced run. Source-file hashes are provenance, not proof of biological validity.

## Status and attribution

Prepared 22 October 2026 with OpenAI/Codex assistance. The sequence analysis and local verification are complete. No biological rescue, clinician review, clinical recommendation or submission is claimed. The earlier proposal remains in the original research workspace and is not overwritten. Author/team attribution and account-specific submission disclosures remain to be completed by the participant; they are not inferred from a computer username.

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

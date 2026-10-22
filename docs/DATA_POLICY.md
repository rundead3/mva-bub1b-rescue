# Repository scope and status

This is a new, isolated repository. The original research workspace was not converted wholesale into Git or uploaded.

Included material:

- A small reference-only snapshot from MANE, RefSeq and UniProt, with source URLs and hashes.
- The accepted Track 1 allele pair carried forward from the completed diagnostic work, with its status and provenance recorded in `data/scenarios.json`.
- Original analysis code, tests, generated sequence consequences and a source-linked research plan.

Excluded material includes raw or aligned reads, VCFs, genotype-scale tables, clinical documents, identifying sample names, local logs, credentials, reference genomes, caches and third-party software archives. No new patient-data upload or recontact is part of this work.

`scripts/check_release.py` uses a conservative allowlist and detects selected raw-file and secret patterns. It is a guard against accidental inclusions, not proof that arbitrary content is deidentified or that Git history is safe. Review the explicitly staged files before publication.

The accepted allele pair originates in the prior restricted-data analysis. Its inclusion does not make the source dataset public. Public reference sequence and diagnostic findings have distinct provenance.

The code is a reference computation, not the full historical variant-discovery pipeline. It takes the accepted Track 1 result as input. Its tests establish software and sequence consistency; functional rescue and clinical safety remain unmeasured.

Hosting is a separate action from local Git initialization. A remote URL should be recorded only after successful creation. No author identity, license grant on behalf of the participant, clinician endorsement or completed submission is inferred. Team attribution, the eventual license and the event's source-data citation remain publication metadata for the participant to finalize.

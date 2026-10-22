# BUB1B sequence feasibility results

Computed from a public reference and the accepted Track 1 BUB1B allele pair. Diagnosis is an input; these results address sequence feasibility, not treatment efficacy.

| Scenario | Native codon | Altered codon | Exact DNA reversion | Exon |
| --- | --- | --- | --- | --- |
| L737_stop | TTA | TGA | G>T | 17 |
| N1002K | AAT | AAG | G>T | 23 |

| Scenario | Chemistry | All compatible nonempty codon edit subsets | Native amino acid restored? |
| --- | --- | --- | --- |
| L737_stop | DNA_ABE_either_strand | CGA (R), CGG (R), TGG (W) | False |
| L737_stop | DNA_CBE_either_strand | TAA (*) | False |
| L737_stop | RNA_ADAR_A_to_I | TGG (W) | False |
| L737_stop | RNA_C_to_U | No compatible base in codon | False |
| N1002K | DNA_ABE_either_strand | AGG (R), GAG (E), GGG (G) | False |
| N1002K | DNA_CBE_either_strand | AAA (K) | False |
| N1002K | RNA_ADAR_A_to_I | AGG (R), GAG (E), GGG (G) | False |
| N1002K | RNA_C_to_U | No compatible base in codon | False |

DNA alphabet is used throughout; RNA T corresponds to U. Products are chemistry-compatible possibilities, not designed or tested editors. Either-strand DNA combinations may require separate targeting events.

The stop has 6 downstream exon junctions. The positional NMD criterion is met; RNA abundance is unmeasured.

Precisely skipping exon 17 removes 141 coding bases. Reading frame preserved: True. Modeled protein: 1003 residues; internal stop: False.

The splice junction splits a codon (phase 1); the modeled fusion codon is GGT (G). Full-sequence comparison verifies a simple deletion of reference residues 715-761 despite the split codon.

Non-allele-selective skipping also creates the shortened protein carrying N1002K on the other copy. Its residual function could be impaired; stop-allele rescue alone would not establish net benefit.

Native Asn1002 can be restored by AAG to AAT (G>T) or AAC (G>C). The latter restores the protein without restoring the original DNA codon. Newer transversion editors are outside the standard ABE/CBE enumeration.

An in-frame shortened protein may fail to fold or function. This result supports investigating the protein product, not prioritizing an exon-skipping medicine.

## Limits

- No PAM, guide, edit-window, off-target, delivery or efficiency assessment.
- Either-strand DNA products are a chemical upper bound and may require separate targeting events.
- No patient RNA, protein, segregation or drug response measured.
- Standard ABE/CBE and two RNA chemistries only; not all emerging editing technologies.
- The accepted Track 1 diagnosis and trans configuration are inputs to this rescue analysis.
- Single-exon skipping models a precise transcript product, not a validated intervention.

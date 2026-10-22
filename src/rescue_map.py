"""Reference-based sequence feasibility, not a clinical efficacy predictor."""
import hashlib
from itertools import product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = dict(zip(map(''.join, product('TCAG', repeat=3)),
                'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
# Effects on the coding strand. DNA editors may act on either strand;
# RNA edits apply to the transcript itself, with I decoded as G.
EDIT_CHEMISTRIES = {
    'DNA_ABE_either_strand': {'A': 'G', 'T': 'C'},
    'DNA_CBE_either_strand': {'C': 'T', 'G': 'A'},
    'RNA_ADAR_A_to_I': {'A': 'G'},
    'RNA_C_to_U': {'C': 'T'},
}


def translate(sequence):
    if len(sequence) % 3 or set(sequence) - set('ACGT'):
        raise ValueError('Expected complete unambiguous codons')
    return ''.join(CODE[sequence[i:i + 3]] for i in range(0, len(sequence), 3))


def edit_products(codon, chemistry):
    """Enumerate all nonempty subsets of compatible edits within one codon.

    This is a chemical upper bound. Either-strand DNA rows may combine edits
    that require separate targeting events, not one guide or editor exposure.
    It does not assess guide, PAM, window,
    bystander edits elsewhere, efficiency, delivery or biological tolerability.
    """
    if len(codon) != 3 or set(codon) - set('ACGT'):
        raise ValueError('Expected a DNA-alphabet codon')
    choices = [(base, chemistry[base]) if base in chemistry else (base,) for base in codon]
    return sorted({''.join(row) for row in product(*choices)} - {codon})


def single_change_native_products(mutant, native_amino_acid):
    rows = []
    for index, original in enumerate(mutant):
        for base in 'ACGT':
            altered = mutant[:index] + base + mutant[index + 1:]
            if base != original and CODE[altered] == native_amino_acid:
                rows.append(dict(codon=altered, codon_position1=index + 1, change=f'{original}>{base}'))
    return rows


def analyze_variant(reference, variant):
    cds = reference['cds']
    p = variant['coding_position1'] - 1
    if not 0 <= p < len(cds) or cds[p] != variant['reference']:
        raise ValueError('Scenario disagrees with public reference')
    begin = p // 3 * 3
    native = cds[begin:begin + 3]
    mutant = native[:p % 3] + variant['alternate'] + native[p % 3 + 1:]
    if (native, mutant) != (variant['expected_reference_codon'], variant['expected_alternate_codon']):
        raise ValueError('Unexpected scenario codons')
    exon = next(x for x in reference['exons'] if x['coding_start1'] <= p + 1 <= x['coding_end1'])
    edits = []
    for name, chemistry in EDIT_CHEMISTRIES.items():
        products = edit_products(mutant, chemistry)
        edits.append(dict(chemistry=name, products=[dict(codon=c, amino_acid=CODE[c]) for c in products],
                          restores_native_codon=native in products,
                          restores_native_amino_acid=any(CODE[c] == CODE[native] for c in products)))
    result = dict(label=variant['label'], coding_position1=p + 1, residue1=p // 3 + 1,
                  reference_codon=native, alternate_codon=mutant,
                  native_amino_acid=CODE[native], altered_amino_acid=CODE[mutant],
                  exact_reversion=f"{variant['alternate']}>{variant['reference']}",
                  exon_number=exon['number'], edit_products=edits,
                  single_base_native_residue_routes=single_change_native_products(mutant, CODE[native]))
    if CODE[mutant] == '*':
        stop_end1 = begin + 3
        junctions = [dict(next_exon=x['number'], next_exon_coding_start1=x['coding_start1'],
                          bases_strictly_between_stop_and_next_exon=x['coding_start1'] - stop_end1 - 1)
                     for x in reference['exons'] if x['coding_start1'] > stop_end1]
        result['downstream_junctions'] = junctions
        result['positional_nmd_criterion'] = any(x['bases_strictly_between_stop_and_next_exon'] > 55 for x in junctions)
        # A precise deletion of one complete internal coding exon is a model,
        # not an assertion that an ASO can achieve that splicing outcome.
        a, b = exon['coding_start1'] - 1, exon['coding_end1']
        if exon['number'] not in (1, len(reference['exons'])):
            shortened = cds[:a] + cds[b:]
            skip = dict(exon=exon['number'], deleted_coding_bases=b - a,
                        reading_frame_preserved=(b - a) % 3 == 0,
                        upstream_codon_phase=a % 3, downstream_codon_phase=b % 3)
            if skip['reading_frame_preserved']:
                protein = translate(shortened)
                # Exon boundaries may split codons; explicitly record the fusion.
                junction_index = a // 3
                skip.update(protein_length=len(protein), internal_stop_present='*' in protein,
                            upstream_last_unchanged_residue1=a // 3,
                            first_retained_original_residue1=b // 3 + 1,
                            junction_codon=shortened[junction_index * 3:junction_index * 3 + 3],
                            junction_residue=protein[junction_index],
                            reference_junction_codon=cds[junction_index * 3:junction_index * 3 + 3],
                            function='unknown; an in-frame product can be nonfunctional')
                first, last = a // 3, b // 3
                if protein == reference['protein'][:first] + reference['protein'][last:]:
                    skip['sequence_verified_simple_deletion'] = dict(first_residue1=first + 1, last_residue1=last,
                        first_reference_amino_acid=reference['protein'][first],
                        last_reference_amino_acid=reference['protein'][last - 1])
            result['single_exon_skip_model'] = skip
    return result


def main():
    paths = [ROOT / 'data/public/bub1b_reference.json', ROOT / 'data/scenarios.json']
    reference, scenarios = [json.loads(p.read_text(encoding='utf-8')) for p in paths]
    if reference['refseq'] != scenarios['transcript'] or translate(reference['cds']) != reference['protein']:
        raise ValueError('Reference identity/translation mismatch')
    if len(reference['cds']) != reference['exons'][-1]['coding_end1']:
        raise ValueError('Reference exon geometry mismatch')
    expected_start = 1
    for exon in reference['exons']:
        if exon['coding_start1'] != expected_start or exon['coding_end1'] < expected_start:
            raise ValueError('Coding intervals must be ordered, contiguous and nonempty')
        expected_start = exon['coding_end1'] + 1
    result = dict(analysis_date='2026-10-22', scope='Sequence feasibility only; no rescue efficacy measured.',
                  diagnostic_status=scenarios['diagnostic_status'],
                  diagnostic_status_basis=scenarios['status_basis'],
                  phase_assumption=scenarios['phase_assumption'],
                  input_sha256={str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  variants=[analyze_variant(reference, v) for v in scenarios['variants']],
                  limits=['No PAM, guide, edit-window, off-target, delivery or efficiency assessment.',
                          'Either-strand DNA products are a chemical upper bound and may require separate targeting events.',
                          'No patient RNA, protein, segregation or drug response measured.',
                          'Standard ABE/CBE and two RNA chemistries only; not all emerging editing technologies.',
                          'The accepted Track 1 diagnosis and trans configuration are inputs to this rescue analysis.',
                          'Single-exon skipping models a precise transcript product, not a validated intervention.'])
    # Non-allele-selective skipping would alter both copies in the supplied trans model.
    stop = next(v for v in result['variants'] if v['altered_amino_acid'] == '*')
    exon = next(x for x in reference['exons'] if x['number'] == stop['exon_number'])
    result['nonselective_skip_products'] = []
    for variant in scenarios['variants']:
        p = variant['coding_position1'] - 1
        sequence = reference['cds'][:p] + variant['alternate'] + reference['cds'][p + 1:]
        skipped = sequence[:exon['coding_start1'] - 1] + sequence[exon['coding_end1']:]
        protein = translate(skipped)
        variant_removed = exon['coding_start1'] <= p + 1 <= exon['coding_end1']
        new_coding_position = None if variant_removed else p + 1 - (exon['coding_end1'] - exon['coding_start1'] + 1 if p + 1 > exon['coding_end1'] else 0)
        result['nonselective_skip_products'].append(dict(allele=variant['label'], protein_length=len(protein),
            scenario_variant_removed=variant_removed,
            retained_variant_original_residue1=None if variant_removed else p // 3 + 1,
            retained_variant_new_residue1=None if variant_removed else (new_coding_position - 1) // 3 + 1,
            retained_variant_amino_acid=None if variant_removed else protein[(new_coding_position - 1) // 3],
            protein_sha256=hashlib.sha256(protein.encode()).hexdigest(),
            internal_stop_present='*' in protein,
            interpretation='Function unknown; skipping the missense-bearing copy could impair residual function.'))
    out = ROOT / 'results'
    out.mkdir(exist_ok=True)
    (out / 'rescue_map.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    lines = ['# BUB1B sequence feasibility results', '',
             'Computed from a public reference and the accepted Track 1 BUB1B allele pair. Diagnosis is an input; these results address sequence feasibility, not treatment efficacy.', '',
             '| Scenario | Native codon | Altered codon | Exact DNA reversion | Exon |',
             '| --- | --- | --- | --- | --- |']
    for v in result['variants']:
        lines.append(f"| {v['label']} | {v['reference_codon']} | {v['alternate_codon']} | {v['exact_reversion']} | {v['exon_number']} |")
    lines += ['', '| Scenario | Chemistry | All compatible nonempty codon edit subsets | Native amino acid restored? |',
              '| --- | --- | --- | --- |']
    for v in result['variants']:
        for row in v['edit_products']:
            products = ', '.join(f"{x['codon']} ({x['amino_acid']})" for x in row['products']) or 'No compatible base in codon'
            lines.append(f"| {v['label']} | {row['chemistry']} | {products} | {row['restores_native_amino_acid']} |")
    skip = stop['single_exon_skip_model']
    lines += ['', 'DNA alphabet is used throughout; RNA T corresponds to U. Products are chemistry-compatible possibilities, not designed or tested editors. Either-strand DNA combinations may require separate targeting events.', '',
              f"The stop has {len(stop['downstream_junctions'])} downstream exon junctions. The positional NMD criterion is met; RNA abundance is unmeasured.", '',
              f"Precisely skipping exon {skip['exon']} removes {skip['deleted_coding_bases']} coding bases. Reading frame preserved: {skip['reading_frame_preserved']}. Modeled protein: {skip['protein_length']} residues; internal stop: {skip['internal_stop_present']}.", '',
              f"The splice junction splits a codon (phase {skip['upstream_codon_phase']}); the modeled fusion codon is {skip['junction_codon']} ({skip['junction_residue']}). Full-sequence comparison verifies a simple deletion of reference residues {skip['sequence_verified_simple_deletion']['first_residue1']}-{skip['sequence_verified_simple_deletion']['last_residue1']} despite the split codon.", '',
              'Non-allele-selective skipping also creates the shortened protein carrying N1002K on the other copy. Its residual function could be impaired; stop-allele rescue alone would not establish net benefit.', '',
              'Native Asn1002 can be restored by AAG to AAT (G>T) or AAC (G>C). The latter restores the protein without restoring the original DNA codon. Newer transversion editors are outside the standard ABE/CBE enumeration.', '',
              'An in-frame shortened protein may fail to fold or function. This result supports investigating the protein product, not prioritizing an exon-skipping medicine.', '',
              '## Limits', ''] + [f'- {x}' for x in result['limits']]
    (out / 'rescue_map.md').write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({v['label']: {'codon': v['alternate_codon'], 'exon': v['exon_number'], 'skip': v.get('single_exon_skip_model')} for v in result['variants']}, indent=2))


if __name__ == '__main__':
    main()

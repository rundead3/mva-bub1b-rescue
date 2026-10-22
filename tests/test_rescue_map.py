import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from rescue_map import analyze_variant, edit_products, translate, single_change_native_products, EDIT_CHEMISTRIES


class RescueMapTests(unittest.TestCase):
    def test_genetic_code_and_invalid_partial_codon(self):
        self.assertEqual(translate('ATGTTATGGTGATAAAAT'), 'MLW**N')
        with self.assertRaises(ValueError):
            translate('AT')

    def test_rna_edit_does_not_apply_reverse_strand_chemistry(self):
        self.assertEqual(edit_products('TGA', EDIT_CHEMISTRIES['RNA_ADAR_A_to_I']), ['TGG'])
        self.assertEqual(edit_products('TGA', EDIT_CHEMISTRIES['RNA_C_to_U']), [])
        self.assertEqual(edit_products('TGA', EDIT_CHEMISTRIES['DNA_CBE_either_strand']), ['TAA'])

    def test_multiple_compatible_edits_are_not_omitted(self):
        self.assertEqual(edit_products('AAG', {'A': 'G'}), ['AGG', 'GAG', 'GGG'])

    def test_native_protein_is_not_limited_to_original_dna_codon(self):
        self.assertEqual({x['codon'] for x in single_change_native_products('AAG', 'N')}, {'AAT', 'AAC'})
        self.assertEqual({x['codon'] for x in single_change_native_products('TGA', 'L')}, {'TTA'})

    def test_reference_and_scenarios(self):
        ref = json.loads((ROOT / 'data/public/bub1b_reference.json').read_text())
        config = json.loads((ROOT / 'data/scenarios.json').read_text())
        self.assertEqual(translate(ref['cds']), ref['protein'])
        self.assertEqual(len(ref['protein']), 1050)
        a, b = [analyze_variant(ref, v) for v in config['variants']]
        self.assertEqual((a['alternate_codon'], b['alternate_codon']), ('TGA', 'AAG'))
        self.assertTrue(all(not x['restores_native_amino_acid'] for v in (a, b) for x in v['edit_products']))
        self.assertEqual(a['downstream_junctions'][-1]['bases_strictly_between_stop_and_next_exon'], 746)
        self.assertEqual(a['single_exon_skip_model']['deleted_coding_bases'], 141)
        self.assertEqual(a['single_exon_skip_model']['protein_length'], 1003)
        deletion = a['single_exon_skip_model']['sequence_verified_simple_deletion']
        self.assertEqual((deletion['first_residue1'], deletion['last_residue1']), (715, 761))
        self.assertEqual(translate(ref['cds'][:2143] + ref['cds'][2284:]), ref['protein'][:714] + ref['protein'][761:])
        mutant = ref['cds'][:3005] + 'G' + ref['cds'][3006:]
        shortened_mutant = translate(mutant[:2143] + mutant[2284:])
        self.assertEqual(shortened_mutant[954], 'K')
        # Refuse a misnumbered/incorrect reference allele instead of silently proceeding.
        invalid = dict(config['variants'][0], reference='A')
        with self.assertRaises(ValueError):
            analyze_variant(ref, invalid)


if __name__ == '__main__':
    unittest.main()

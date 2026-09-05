"""Regression evidence for the legacy audit, including intentional old defects."""
import importlib.util
from fractions import Fraction as F
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location(
    "version1_reproduction", Path(__file__).resolve().parents[1] / "experiments" / "reproduce_version1.py")
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)


class Version1AuditTests(unittest.TestCase):
    def test_sickle_screenshot_and_missing_mass_are_preserved(self):
        history = legacy.legacy_sickle("AS", "AS", 3)
        self.assertEqual(history, [[F(1,4),F(1,2),F(1,4)],
                                   [F(9,16),F(3,8),F(1,16)],
                                   [F(3,4),F(3,16),F(0)]])
        self.assertEqual(sum(history[-1]), F(15,16))

    def test_sickle_screenshot_is_not_plain_matrix_power(self):
        initial = legacy.punnett("AS", "AS", legacy.SICKLE)
        plain = legacy.matvec(legacy.M_SICKLE, legacy.matvec(legacy.M_SICKLE, initial))
        self.assertEqual(plain, [F(49,64),F(7,32),F(1,64)])
        self.assertNotEqual(plain, legacy.legacy_sickle("AS", "AS", 3)[-1])

    def test_abo_examples_use_all_pairings(self):
        aa_oo = legacy.legacy_random_mating("AA","OO",3)
        self.assertEqual(aa_oo[-1], [F(1,4),F(1,2),F(0),F(0),F(0),F(1,4)])
        ao_bo = legacy.legacy_random_mating("AO","BO",4)
        self.assertEqual(ao_bo[-1], [F(1,16),F(1,4),F(1,16),F(1,4),F(1,8),F(1,4)])
        self.assertTrue(all(sum(row)==1 for row in aa_oo+ao_bo))

    def test_abo_rh_p49_stationary_from_generation_one(self):
        history = legacy.legacy_random_mating("AO/Dd","AO/Dd",3,joint=True)
        self.assertTrue(all(x == history[0] for x in history))
        self.assertEqual(legacy.phenotypes(history[-1],True),
                         [F(9,16),F(3,16),F(0),F(0),F(0),F(0),F(3,16),F(1,16)])

    def test_abo_rh_p50_independent_allele_truth(self):
        # AA/dd x BO/Dd yields A=.5, B=.25, O=.25 and d=.75.
        truth_abo = [F(1,4),F(1,4),F(1,16),F(1,8),F(1,4),F(1,16)]
        truth_rh = [F(1,16),F(3,8),F(9,16)]
        final = legacy.legacy_random_mating("AA/dd","BO/Dd",7,joint=True)[-1]
        self.assertEqual(final, [a*r for a in truth_abo for r in truth_rh])

    def test_displayed_joint_matrix_is_self_cross_not_ranked_pairs(self):
        for c,g in enumerate(legacy.JOINT):
            self.assertEqual([row[c] for row in legacy.M_JOINT_DISPLAYED], legacy.punnett18(g,g))
        different = sum(legacy.M_JOINT_DISPLAYED[r][c] != legacy.punnett18(a,b)[r]
                        for c,(a,b,_) in enumerate(legacy.PAPER_TOP18) for r in range(18))
        self.assertEqual(different,114)

    def test_complete_pair_count_and_mass(self):
        pairs = legacy.pair_weights(legacy.ABO,[F(1,6)]*6)
        self.assertEqual(len(pairs),21)
        self.assertEqual(sum(w for _,_,w in pairs),1)
        self.assertEqual(len(legacy.pair_weights(legacy.JOINT,[F(1,18)]*18)),171)

    def test_coverage_rounding_provenance(self):
        audit,_ = legacy.coverage_audit()
        exact = audit["variants"]["paper_alleles_exact"]
        raw = audit["variants"]["paper_p32_pair_formula_inputs_raw"]
        self.assertAlmostEqual(exact["top6_percent_raw_mass"],83.4764828157,places=10)
        self.assertAlmostEqual(exact["top18_joint_percent_raw_mass"],67.83473349212181,places=11)
        self.assertAlmostEqual(raw["top6_percent_raw_mass"],83.6205,places=10)
        self.assertAlmostEqual(raw["total_unordered_pair_mass"],1.002001,places=12)
        self.assertEqual(audit["paper_top18_printed_row_sum_percent"],F("67.84"))


if __name__ == "__main__":
    unittest.main()

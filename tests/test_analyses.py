#! /usr/bin/env python

import os
import math
import random
import unittest
import json
import subprocess
import dendropy
from dendropy.utility import processio

if __name__ == "__main__":
    # For:
    #   python tests/test_probabilities.py
    import _pathmap
else:
    # For:
    #   python setup.py test
    #   python -m unittest discover
    #   python -m unittest tests/test_probabilities.py
    from . import _pathmap

class ConstrainedPartitionsTestCase(unittest.TestCase):

    def execute_analysis(self,
            config_path,
            tree_path):
        cmd = [
                os.path.join(_pathmap.BIN_DIR, "delineate-estimate-species-partition.py"),
                "-c", config_path,
                "-t", tree_path,
                "-I",
              ]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,)
        stdout, stderr = processio.communicate(p)
        return json.loads(stdout)

    def _load_partitions(self, d):
        partitions = {}
        for part_d in d["partitions"]:
            species_leafsets = frozenset(frozenset(s) for s in part_d["species_leafsets"])
            partitions[species_leafsets] = {}
            for val_key in (
                    "log_probability",
                    "probability",
                    "probability_given_constraints",
                    "log_probability_given_constraints",
                    "is_in_confidence_interval",
                    "cumulative_probability",
                    "cumulative_probability_given_constraints",
                    ):
                partitions[species_leafsets][val_key] = part_d[val_key]
        return partitions

    def _check_analysis(self,
            config_path,
            tree_path,
            expected_results_path):
        with open(expected_results_path, "r") as src:
            expected_results = json.load(src)
        observed_results = self.execute_analysis(
            config_path=config_path,
            tree_path=tree_path)
        expected_partitions = self._load_partitions(expected_results)
        observed_partitions = self._load_partitions(observed_results)
        for part_key in expected_partitions:
            exp_vals = expected_partitions[part_key]
            obs_vals = observed_partitions[part_key]
            for val_key in exp_vals:
                exp_val = exp_vals[val_key]
                obs_val = obs_vals[val_key]
                if type(exp_val) == float:
                    self.assertAlmostEqual(exp_val, obs_val, 8)
                else:
                    self.assertEqual(exp_val, obs_val)

        for key in [
                "speciation_completion_rate",
                "speciation_completion_rate_estimate_lnl",
                "num_partitions",
                "max_log_likelihood",
                "num_partitions_in_confidence_interval",
                ]:
            exp_val = expected_results[key]
            obs_val = observed_results[key]
            if type(exp_val) == float:
                self.assertAlmostEqual(exp_val, obs_val, 8)
            else:
                self.assertEqual(exp_val, obs_val)

    def test_constrained_partitions_small(self):
        test_file_dir = os.path.join(_pathmap.TESTS_DATA_DIR, "constrained-partitions-small/s1-master-aa56716")
        test_filename_stems = [
                "run1_spr0.001_.0001",
                "run1_spr0.005_.0001",
                "run1_spr0.010_.0001",
                "run1_spr0.050_.0001",
                "run1_spr0.010_.0001",
                ]
        for test_filename_stem in test_filename_stems:
            config_path = os.path.join(test_file_dir, test_filename_stem + ".json")
            tree_path = os.path.join(test_file_dir, test_filename_stem + ".nex")
            expected_results_path = os.path.join(test_file_dir, test_filename_stem + ".partition-probs.json")
            self._check_analysis(
                    config_path=config_path,
                    tree_path=tree_path,
                    expected_results_path=expected_results_path)

if __name__ == "__main__":
    unittest.main()










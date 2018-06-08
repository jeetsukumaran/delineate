#! /usr/bin/env python

import os
import math
import random
import unittest
import json
import dendropy
from delineate import model
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

class LineageTreeBasicMarginalSpeciesProbabilities(unittest.TestCase):

    # Expected results have been independently calculated by JS and MTH, manually
    # and by script.

    def setUp(self):
        s = "[&R]((a:0.20,b:0.20)i7:0.70,(c:0.55,(d:0.30,e:0.30)i6:0.25)i8:0.35);"
        self.tree = model.LineageTree.get(
                data=s,
                schema="newick",
                )
        node_edge_labels = {
                'a': 'v1',
                'b': 'v2',
                'c': 'v3',
                'd': 'v4',
                'e': 'v5',
                'i6': 'v6',
                'i7': 'v7',
                'i8': 'v8',
                None: None,
                }
        for node in self.tree:
            if node.taxon is not None:
                node.edge.label = node_edge_labels[node.taxon.label]
            elif node.label is not None:
                node.edge.label = node_edge_labels[node.label]
        self.monophyletic_multitaxon_clade_test_cases = {
                "ab": 0.020668831269136167,
                "cde": 0.02022241120317388,
                "de": 0.005151702704502191,
        }
        self.single_taxon_clade_test_cases = {
                "a": 0.00407485155242,
                "b": 0.00407485155242,
                "c": 0.0110430425835,
                "d": 0.006013039105285815,
                "e": 0.006013039105285815,
        }
        self.non_root_spanning_nonmonophyletic_multitaxon_clade_test_cases = {
                "ce": 0.00012169919972082542,
        }
        self.root_spanning_nonmonophyletic_multitaxon_clade_test_cases = {
                "ace": 2.2783942169371514e-05,
                "abc": 0.004826167453992684,
        }
        self.tree.speciation_completion_rate = 0.02

    def testValidMonophyleticMultitaxonCladeProbabilities(self):
        for tax_labels in self.monophyletic_multitaxon_clade_test_cases:
            self.assertAlmostEqual(
                    self.tree.calc_marginal_probability_of_species(tax_labels),
                    self.monophyletic_multitaxon_clade_test_cases[tax_labels], 8)

    def testValidSingleTaxonCladeProbabilities(self):
        for tax_labels in self.single_taxon_clade_test_cases:
            self.assertAlmostEqual(
                    self.tree.calc_marginal_probability_of_species(tax_labels),
                    self.single_taxon_clade_test_cases[tax_labels], 8)

    def testValidRootSpanningNonMonophyleticMultitaxonCladeProbabilities(self):
        for tax_labels in self.root_spanning_nonmonophyletic_multitaxon_clade_test_cases:
            self.assertAlmostEqual(self.tree.calc_marginal_probability_of_species(tax_labels),
                    self.root_spanning_nonmonophyletic_multitaxon_clade_test_cases[tax_labels], 8)

    def testValidNonRootSpanningNonMonophyleticMultitaxonCladeProbabilities(self):
        for tax_labels in self.non_root_spanning_nonmonophyletic_multitaxon_clade_test_cases:
            self.assertAlmostEqual(self.tree.calc_marginal_probability_of_species(tax_labels),
                    self.non_root_spanning_nonmonophyletic_multitaxon_clade_test_cases[tax_labels], 8)

class LineageTreeMultiMarginalSpeciesProbabilities(unittest.TestCase):

    # Tests to see if changing the branch lengths and/or speciation rate will
    # result in correct probs
    def test_probs(self):
        with open(os.path.join(_pathmap.TESTS_DATA_DIR, "marginal_probability_of_species.json")) as src:
            test_ref = json.load(src)
        for test_tree_set in test_ref:
            taxon_namespace = dendropy.TaxonNamespace(test_tree_set["taxon_namespace"])
            tree = model.LineageTree.get(
                    data=test_tree_set["tree"],
                    schema="newick",
                    taxon_namespace=taxon_namespace,
                    )
            tree.encode_bipartitions()
            for brlen_config in test_tree_set["branch_length_configurations"]:
                for split_bitmask, br_len in brlen_config["branch_lengths"].items():
                    split_bitmask = int(split_bitmask)
                    assert split_bitmask in tree.split_bitmask_edge_map, split_bitmask
                    tree.split_bitmask_edge_map[split_bitmask].length = br_len
                for speciation_rate_config in brlen_config["speciation_rate_configurations"]:
                    tree.speciation_completion_rate = speciation_rate_config["speciation_rate"]
                    for species_configuration in speciation_rate_config["species_configurations"]:
                        species_labels = species_configuration["species"]
                        expected_probability = species_configuration["probability"]
                        obs_probability = tree.calc_marginal_probability_of_species(species_labels)
                        self.assertAlmostEqual(expected_probability, obs_probability, 8)

class LineageTreeJointSpeciesProbabilities(unittest.TestCase):

    # Tests to see if changing the branch lengths and/or speciation rate will
    # result in correct probs
    def test_probs(self):
        with open(os.path.join(_pathmap.TESTS_DATA_DIR, "joint_probability_of_species.json")) as src:
            test_ref = json.load(src)
        for test_tree_set in test_ref:
            taxon_namespace = dendropy.TaxonNamespace(test_tree_set["taxon_namespace"])
            tree = model.LineageTree.get(
                    data=test_tree_set["tree"],
                    schema="newick",
                    taxon_namespace=taxon_namespace,
                    )
            tree.encode_bipartitions()
            for brlen_config in test_tree_set["branch_length_configurations"]:
                for split_bitmask, br_len in brlen_config["branch_lengths"].items():
                    split_bitmask = int(split_bitmask)
                    assert split_bitmask in tree.split_bitmask_edge_map, split_bitmask
                    tree.split_bitmask_edge_map[split_bitmask].length = br_len
                for speciation_rate_config in brlen_config["speciation_rate_configurations"]:
                    tree.speciation_completion_rate = speciation_rate_config["speciation_rate"]
                    for spconf_idx, species_configuration in enumerate(speciation_rate_config["species_configurations"]):
                        species_leafset_labels = frozenset([frozenset(s) for s in species_configuration["species"]])
                        expected_probability=species_configuration["probability"]
                        obs_probability=tree.calc_joint_probability_of_species(species_leafset_labels)
                        self.assertAlmostEqual(expected_probability, obs_probability, 8)

class LineageTreeJointSpeciesProbabilitiesWithConstraints(unittest.TestCase):

    def test_probs_with_fully_specified_constraints_on_new_tree(self):
        with open(os.path.join(_pathmap.TESTS_DATA_DIR, "joint_probability_of_species.json")) as src:
            test_ref = json.load(src)
        for test_tree_set in test_ref:
            taxon_namespace = dendropy.TaxonNamespace(test_tree_set["taxon_namespace"])
            for brlen_config in test_tree_set["branch_length_configurations"]:
                for speciation_rate_config in brlen_config["speciation_rate_configurations"]:
                    for spconf_idx, species_configuration in enumerate(speciation_rate_config["species_configurations"]):
                        species_leafset_labels=frozenset([frozenset(s) for s in species_configuration["species"]])
                        tree = model.LineageTree.get(
                                data=test_tree_set["tree"],
                                schema="newick",
                                taxon_namespace=taxon_namespace,
                                )
                        tree.encode_bipartitions()
                        for split_bitmask, br_len in brlen_config["branch_lengths"].items():
                            split_bitmask = int(split_bitmask)
                            assert split_bitmask in tree.split_bitmask_edge_map, split_bitmask
                            tree.split_bitmask_edge_map[split_bitmask].length = br_len
                        tree.speciation_completion_rate = speciation_rate_config["speciation_rate"]
                        tree.set_node_constraints(species_leafset_labels=species_leafset_labels)
                        expected_probability=species_configuration["probability"]
                        obs_probability=tree.calc_joint_probability_of_species(species_leafset_labels)
                        self.assertAlmostEqual(expected_probability, obs_probability, 8)

#     def test_probs_with_constraints_with_same_tree(self):
#         with open(os.path.join(_pathmap.TESTS_DATA_DIR, "joint_probability_of_species.json")) as src:
#             test_ref = json.load(src)
#         for test_tree_set in test_ref:
#             taxon_namespace = dendropy.TaxonNamespace(test_tree_set["taxon_namespace"])
#             tree = model.LineageTree.get(
#                     data=test_tree_set["tree"],
#                     schema="newick",
#                     taxon_namespace=taxon_namespace,
#                     )
#             tree.encode_bipartitions()
#             for brlen_config in test_tree_set["branch_length_configurations"]:
#                 for split_bitmask, br_len in brlen_config["branch_lengths"].items():
#                     split_bitmask = int(split_bitmask)
#                     assert split_bitmask in tree.split_bitmask_edge_map, split_bitmask
#                     tree.split_bitmask_edge_map[split_bitmask].length = br_len
#                 for speciation_rate_config in brlen_config["speciation_rate_configurations"]:
#                     tree.speciation_completion_rate = speciation_rate_config["speciation_rate"]
#                     for spconf_idx, species_configuration in enumerate(speciation_rate_config["species_configurations"]):
#                         species_leafset_labels=frozenset([frozenset(s) for s in species_configuration["species"]])
#                         expected_probability=species_configuration["probability"]
#                         obs_probability=tree.calc_joint_probability_of_species(species_leafset_labels)
#                         self.assertAlmostEqual(expected_probability, obs_probability, 8)


if __name__ == "__main__":
    unittest.main()




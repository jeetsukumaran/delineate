#! /usr/bin/env python

import math
import random
import unittest
from delineate import model

class LineageTreeSpeciationProbabilities(unittest.TestCase):

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
        rng = random.Random()
        species_conversion_rate = model.GammaDistributedParameter(rng=rng, shape=0.1, scale=0.01, initial_value=0.02)
        self.tree.build(rng=rng, species_conversion_rate=species_conversion_rate)
        self.expected_probs = {
            'v1': (0.996007989344, 0.00399201065601),
            'v2': (0.996007989344, 0.00399201065601),
            'v3': (0.989060278775, 0.0109397212246),
            'v4': (0.994017964054, 0.00598203594606),
            'v5': (0.994017964054, 0.00598203594606),
            'v6': (0.995012479193, 0.00498752080732),
            'v7': (0.986097544263, 0.0139024557371),
            'v8': (0.993024442933, 0.00697555706676),
        }

    def testEdgeProbs(self):
        for nd in self.tree:
            if nd.edge.label is None:
                continue
            exp_probability_of_no_speciation, exp_probability_of_any_speciation = self.expected_probs[nd.edge.label]
            self.assertAlmostEqual(exp_probability_of_no_speciation, nd.edge.probability_of_no_speciation, 8)
            self.assertAlmostEqual(exp_probability_of_any_speciation, nd.edge.probability_of_any_speciation, 8)

    def testEdgeProbsNoSpeciation(self):
        for nd in self.tree:
            if nd.edge.label is None:
                continue
            exp_probability_of_no_speciation, exp_probability_of_any_speciation = self.expected_probs[nd.edge.label]
            self.assertAlmostEqual(exp_probability_of_no_speciation, nd.edge.probability_of_no_speciation, 8)

    def testEdgeProbsAnySpeciation(self):
        for nd in self.tree:
            if nd.edge.label is None:
                continue
            exp_probability_of_no_speciation, exp_probability_of_any_speciation = self.expected_probs[nd.edge.label]
            self.assertAlmostEqual(exp_probability_of_any_speciation, nd.edge.probability_of_any_speciation, 8)

    def testEdgeLogProbs(self):
        for nd in self.tree:
            if nd.edge.label is None:
                continue
            exp_log_probability_of_no_speciation, exp_log_probability_of_any_speciation = [math.log(i) for i in self.expected_probs[nd.edge.label]]
            self.assertAlmostEqual(exp_log_probability_of_no_speciation, nd.edge.log_probability_of_no_speciation, 8)
            self.assertAlmostEqual(exp_log_probability_of_any_speciation, nd.edge.log_probability_of_any_speciation, 8)

    def testEdgeLogProbsNoSpeciation(self):
        for nd in self.tree:
            if nd.edge.label is None:
                continue
            exp_log_probability_of_no_speciation, exp_log_probability_of_any_speciation = [math.log(i) for i in self.expected_probs[nd.edge.label]]
            self.assertAlmostEqual(exp_log_probability_of_any_speciation, nd.edge.log_probability_of_any_speciation, 8)

    def testEdgeLogProbsAnySpeciation(self):
        for nd in self.tree:
            if nd.edge.label is None:
                continue
            exp_log_probability_of_no_speciation, exp_log_probability_of_any_speciation = [math.log(i) for i in self.expected_probs[nd.edge.label]]
            self.assertAlmostEqual(exp_log_probability_of_any_speciation, nd.edge.log_probability_of_any_speciation, 8)

    def testEdgeLogProbsAndProbs(self):
        for nd in self.tree:
            if nd.edge.label is None:
                continue
            exp_probability_of_no_speciation, exp_probability_of_any_speciation = self.expected_probs[nd.edge.label]
            exp_log_probability_of_no_speciation, exp_log_probability_of_any_speciation = [math.log(i) for i in self.expected_probs[nd.edge.label]]
            self.assertAlmostEqual(exp_log_probability_of_no_speciation, nd.edge.log_probability_of_no_speciation, 8)
            self.assertAlmostEqual(exp_log_probability_of_any_speciation, nd.edge.log_probability_of_any_speciation, 8)
            self.assertAlmostEqual(exp_probability_of_no_speciation, nd.edge.probability_of_no_speciation, 8)
            self.assertAlmostEqual(exp_probability_of_any_speciation, nd.edge.probability_of_any_speciation, 8)

    def testEdgeProbsAndLogProbs(self):
        for nd in self.tree:
            if nd.edge.label is None:
                continue
            exp_probability_of_no_speciation, exp_probability_of_any_speciation = self.expected_probs[nd.edge.label]
            exp_log_probability_of_no_speciation, exp_log_probability_of_any_speciation = [math.log(i) for i in self.expected_probs[nd.edge.label]]
            self.assertAlmostEqual(exp_probability_of_no_speciation, nd.edge.probability_of_no_speciation, 8)
            self.assertAlmostEqual(exp_probability_of_any_speciation, nd.edge.probability_of_any_speciation, 8)
            self.assertAlmostEqual(exp_log_probability_of_no_speciation, nd.edge.log_probability_of_no_speciation, 8)
            self.assertAlmostEqual(exp_log_probability_of_any_speciation, nd.edge.log_probability_of_any_speciation, 8)

class LineageTreeSpeciesProbabilities(unittest.TestCase):

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
        rng = random.Random()
        species_conversion_rate = model.GammaDistributedParameter(rng=rng, shape=0.1, scale=0.01, initial_value=0.02)
        self.tree.build(rng=rng, species_conversion_rate=species_conversion_rate)
        self.monophyletic_multitaxon_clade_test_cases = {
                "ab": 0.020668831269136167,
                "cde": 0.02022241120317388,
                "de": 0.005151702704502191,
        }
        self.single_taxon_clade_test_cases = {
                "a": 0.00399201065601,
                "b": 0.00399201065601,
                "c": 0.0109397212246,
                "d": 0.00598203594606,
                "e": 0.00598203594606,
        }

    def testValidMonophyleticMultitaxonCladeProbabilities(self):
        for tax_labels in self.monophyletic_multitaxon_clade_test_cases:
            taxa = self.tree.taxon_namespace.get_taxa(labels=tax_labels)
            assert len(taxa) == len(tax_labels)
            self.assertAlmostEqual(self.tree.probability_of_monophyletic_multitaxon_clade_good_species(taxa), self.monophyletic_multitaxon_clade_test_cases[tax_labels], 8)

    def testValidSingleTaxonCladeProbabilities(self):
        for tax_label in self.single_taxon_clade_test_cases:
            taxon = self.tree.taxon_namespace.get_taxon(label=tax_label)
            assert taxon is not None
            self.assertAlmostEqual(self.tree.probability_of_single_taxon_good_species(taxon), self.single_taxon_clade_test_cases[tax_label], 8)

if __name__ == "__main__":
    unittest.main()




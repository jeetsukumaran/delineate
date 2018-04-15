#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generate data for tests: for each distinct topology, generates multiple
sub-replicates with varying branch lengths and speciation rates. This will
allow the same Tree object to be re-scored in tests to ensure no artifacts
previous calculation arise.
"""

import sys
import os
import random
import subprocess
import collections
import json
import dendropy
from dendropy.simulate import treesim

script_path = os.path.dirname(__file__)

def randomize_brlens(tree, rng):
    """
    Randomize branch lengths on a tree while maintaining ultrametricity.
    """
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            nd.age = 0.0
        else:
            max_child_age = max(ch.age for ch in nd.child_node_iter())
            nd.age = max_child_age + rng.uniform(0,1)
    tree.set_edge_lengths_from_node_ages()
    tree.seed_node.edge.length = rng.uniform(0,1)

def main():
    """
    Main CLI handler.
    """
    # num_tax = 6
    # num_trees = 3
    # num_branch_length_variants = 3
    # speciation_rates = [0.01, 0.05, 0.2]
    num_tax = 3
    num_trees = 1
    num_branch_length_variants = 1
    speciation_rates = [0.01]
    rng = random.Random()
    working_filepath = ".temp-test-data-tree"
    assert num_tax <= 26
    joint_probability_test_data = []
    marginal_probability_test_data = []
    for tree_idx in range(num_trees):
        taxon_namespace = dendropy.TaxonNamespace(chr(i+97) for i in range(num_tax))
        tree = treesim.birth_death_tree(
                birth_rate=0.02,
                death_rate=0.0,
                num_extant_tips=num_tax,
                taxon_namespace=taxon_namespace)
        tree.encode_bipartitions()
        main_entries = {
                "joint": collections.OrderedDict(),
                "marginal": collections.OrderedDict(),
        }
        tree_string = tree.as_string("newick").replace("\n", "")
        for main_entry in main_entries.values():
            main_entry["tree"] = tree_string
            main_entry["branch_length_configurations"] = []
        assert len(taxon_namespace) == num_tax
        leaf_count = 0
        taxa = [taxon for taxon in taxon_namespace]
        rng.shuffle(taxa)
        leaves = [nd for nd in tree.leaf_node_iter()]
        rng.shuffle(leaves)
        assert len(leaves) == len(taxa)
        for nd, taxon in zip(leaves, taxa):
            nd.taxon = taxon
        for brlen_variant_idx in range(num_branch_length_variants):
            randomize_brlens(tree, rng)
            with open(working_filepath, "w") as dest:
                tree.write(file=dest, schema="newick")
                dest.flush()
            for main_entry in main_entries.values():
                main_entry["branch_length_configurations"].append(collections.OrderedDict())
                main_entry["branch_length_configurations"][-1]["branch_lengths"] = collections.OrderedDict([ (edge.split_bitmask, edge.length) for edge in tree.preorder_edge_iter() ])
                main_entry["branch_length_configurations"][-1]["speciation_rate_configurations"] = []
            for speciation_rate in speciation_rates:
                for main_entry in main_entries.values():
                    main_entry["branch_length_configurations"][-1]["speciation_rate_configurations"].append(collections.OrderedDict())
                    main_entry["branch_length_configurations"][-1]["speciation_rate_configurations"][-1]["speciation_rate"] = speciation_rate
                    main_entry["branch_length_configurations"][-1]["speciation_rate_configurations"][-1]["species_probabilities"] = []
                cmd = [os.path.abspath(os.path.join(script_path, "check.sh")),
                        working_filepath,
                        str(speciation_rate)]
                # Python 3.6 or higher should just use 'subprocess.run()'
                #   subprocess.run(args, *, stdin=None, input=None, stdout=None, stderr=None, shell=False, cwd=None, timeout=None, check=False, encoding=None, errors=None)
                # p = subprocess.Popen(cmd,
                #         stdout=subprocess.PIPE,
                #         stdin=subprocess.PIPE,
                #         )
                # stdout, stderr = p.communicate()
                # if p.returncode:
                #     sys.exit("{} failures reported".format(p.returncode))
                p = subprocess.run(cmd,
                        stdout=subprocess.PIPE,
                        stderr=None,
                        universal_newlines=True,
                        # encoding="utf-8", # not needed if 'universal_newlines' specified?
                        )
                if p.returncode:
                    sys.exit("{} failures reported".format(p.returncode))
                for row in p.stdout.split("\n"):
                    if not row:
                        continue
                    cols = row.split("\t")
                    subentry = collections.OrderedDict((
                            # ("speciation_rate", speciation_rate),
                            ("species", None),
                            ("probability_type", cols[2]),
                            ("probability_value", float(cols[4])),
                            ))
                    if cols[2] == "joint":
                        subentry["species"] = [sp.split(",") for sp in cols[3].split(";")]
                        main_entries["joint"]["branch_length_configurations"][-1]["speciation_rate_configurations"][-1]["species_probabilities"].append(subentry)
                    elif cols[2] == "marginal":
                        subentry["species"] = cols[3].split(";")
                        main_entries["marginal"]["branch_length_configurations"][-1]["speciation_rate_configurations"][-1]["species_probabilities"].append(subentry)
                    else:
                        raise ValueError(cols[2])
        marginal_probability_test_data.append(main_entries["marginal"])
        joint_probability_test_data.append(main_entries["joint"])
    with open(os.path.join("out", "marginal_probability.json"), "w") as dest:
        json.dump(marginal_probability_test_data, dest, indent=4, separators=(',', ': '))
    with open(os.path.join("out", "joint_probability.json"), "w") as dest:
        json.dump(joint_probability_test_data, dest, indent=4, separators=(',', ': '))

if __name__ == "__main__":
    main()

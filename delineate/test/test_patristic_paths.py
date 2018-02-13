#! /usr/bin/env python

import unittest
from delineate import model

class PatristicPathTest(unittest.TestCase):

    def test1(self):
        s = "[&R]((a,b)i7,(c,(d,e)i6)i8);"
        tree = model.LineageTree.get(data=s, schema="newick")
        node_edge_labels = {
                'a': 'v1',
                'b': 'v2',
                'c': 'v3',
                'd': 'v4',
                'e': 'v5',
                'i7': 'v7',
                'i6': 'v6',
                'i8': 'v8',
                }
        for node in tree:
            if node.taxon is not None:
                node.edge.label = node_edge_labels[node.taxon.label]
            elif node.label is not None:
                node.edge.label = node_edge_labels[node.label]
        pdm = tree.phylogenetic_distance_matrix(is_store_path_edges=True)
        expected = {
                ('a', 'a'):  [],
                ('a', 'b'):  ['v1', 'v2'],
                ('a', 'c'):  ['v1', 'v7', 'v8', 'v3'],
                ('a', 'd'):  ['v1', 'v7', 'v8', 'v6', 'v4'],
                ('a', 'e'):  ['v1', 'v7', 'v8', 'v6', 'v5'],
                ('b', 'a'):  ['v2', 'v1'],
                ('b', 'b'):  [],
                ('b', 'c'):  ['v2', 'v7', 'v8', 'v3'],
                ('b', 'd'):  ['v2', 'v7', 'v8', 'v6', 'v4'],
                ('b', 'e'):  ['v2', 'v7', 'v8', 'v6', 'v5'],
                ('c', 'a'):  ['v3', 'v8', 'v7', 'v1'],
                ('c', 'b'):  ['v3', 'v8', 'v7', 'v2'],
                ('c', 'c'):  [],
                ('c', 'd'):  ['v3', 'v6', 'v4'],
                ('c', 'e'):  ['v3', 'v6', 'v5'],
                ('d', 'a'):  ['v4', 'v6', 'v8', 'v7', 'v1'],
                ('d', 'b'):  ['v4', 'v6', 'v8', 'v7', 'v2'],
                ('d', 'c'):  ['v4', 'v6', 'v3'],
                ('d', 'd'):  [],
                ('d', 'e'):  ['v4', 'v5'],
                ('e', 'a'):  ['v5', 'v6', 'v8', 'v7', 'v1'],
                ('e', 'b'):  ['v5', 'v6', 'v8', 'v7', 'v2'],
                ('e', 'c'):  ['v5', 'v6', 'v3'],
                ('e', 'd'):  ['v5', 'v4'],
                ('e', 'e'):  [],
        }
        for t1 in tree.taxon_namespace:
            for t2 in tree.taxon_namespace:
                obs_edges1 = pdm.path_edges(t1, t2)
                obs_edges1_labels = [e.label for e in obs_edges1]
                # print("{}, {}: {}".format(t1.label, t2.label, obs_edges1_labels))
                self.assertEqual(expected[(t1.label, t2.label)], obs_edges1_labels)

if __name__ == "__main__":
    unittest.main()




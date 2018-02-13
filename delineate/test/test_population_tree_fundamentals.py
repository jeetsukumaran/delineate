#! /usr/bin/env python

import unittest
from delineate import model

class PhylogeneticPathTest(unittest.TestCase):

    def test1(self):
        s = "[&R]((a,b)i7,(c,(d,e)i6)i8);"
        tree = model.LineageTree.get(data=s, schema="newick")
        self.assertIs(type(tree), model.LineageTree)
        for nd in tree:
            self.assertIs(type(nd), model.LineageNode)
            self.assertIs(type(nd.edge), model.LineageEdge)

if __name__ == "__main__":
    unittest.main()



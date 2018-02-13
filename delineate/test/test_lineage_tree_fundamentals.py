#! /usr/bin/env python

import random
import unittest
from delineate import model

class LineageTreeFundamentals(unittest.TestCase):

    def test1(self):
        s = "[&R]((a,b)i7,(c,(d,e)i6)i8);"
        rng = random.Random()
        shape = rng.uniform(1e-6, 10)
        scale = rng.uniform(1e-6, 10)
        species_conversion_rate = model.GammaDistributedParameter(rng=rng, shape=shape, scale=scale)
        tree = model.LineageTree.get(
                data=s,
                schema="newick",
                )
        tree.parameterize(rng=rng, species_conversion_rate=species_conversion_rate)
        self.assertIs(type(tree), model.LineageTree)
        self.assertIs(tree.rng, rng)
        self.assertIs(tree.species_conversion_rate, species_conversion_rate)
        for nd in tree:
            self.assertIs(type(nd), model.LineageNode)
            self.assertIs(nd.tree, tree)
            self.assertIs(type(nd.edge), model.LineageEdge)
            self.assertIs(nd.edge.tree, tree)

if __name__ == "__main__":
    unittest.main()



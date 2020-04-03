

1.  A preliminary StarBeast2 analysis, treating all geographic localities as distinct "species" (really, populations).
    The purpose of this is to provide a guide tree for the BPP A10 analysis, which we will use to identify the true population lineages.

2.  BPP analysis, treating all geographic localities as distinct "species" (really, populations).
    This analysis will collapse all lineages between which no restriction in gene flow can be detected under the multipopulation coalescent (or Multispecies Coalescent, MSC) model.
    The collapsed lineages are our actual Wright-Fisher populations that shall form the fundamental units of the population lienage tree.

    We carry out this analysis separate for each (actual) nominal species, as BPP cannot handle the full set of lineages at once.
    Within each nominal species in the original dataset, there are several geographic localities identified.
    We assume before the analysis that it is possible that each of these geographic localities constitute an independent Wright-Fisher population.
    We will carry out a BPP A10 analysis, i.e., a "species" delimitation analysis, but what we are trying to do here is to actually delimit *populations*.
    We will use BPP to

3.  StarBeast2 analysis, with the population lineages identified in the BPP step before functioning as "species".
    This analysis will give us the population tree (i.e., the "species tree", or the MSC ultrametric tree), along with the crucial (relative) branch lengths, for the DELINEATE analysis.

4.  A set of DELINEATE analyses, taking the population tree from the previous step as input.

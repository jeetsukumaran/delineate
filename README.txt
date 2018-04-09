
## Examples

As brute-force calculation of the joint species scenarios by enumerating
all bipartitions of branches into the categories of:
  1. "had 0 speciation events", and 
  2. "had at least one speciation event"
one can run:

    python delineate/joint_from_all_scenarios.py example/five_leaf.tre 0.02

where the 0.02 is a scalar that is multiplied by each branch length to get
    the expected number of good speciation events along that branch.
This script will generat the full joint distribution, and then summarize
    it by producing the marginal probabilities for any species delimitation.


The more efficient, dynamic programming approach for getting the joint probabilities
is found in `efficient_joint_prob.py`:

    python delineate/efficient_joint_prob.py example/five_leaf.tre 0.02


If you just want the probability of one subset of leaves being recognized
    as a single species (distinct from all other tips) then you can run:

    python delineate/marginal_prob.py example/five_leaf.tre .02 a c e

to get the probability that `{a, c, e}` is a species.


## Tests

    ./check.sh example/five_leaf.tre 0.02


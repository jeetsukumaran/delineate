####################
Fundamental Concepts
####################

"Partitions"
============

The primary output of |delineate| is essentially calculation of the probabilities of different *partitions* of the population lineages in the study into distinct species.
Each partition represents a different collection of assignments of species identities to each of the population lineages.
The term "partition" here comes from basic set theory: it is a division of a set of elements into a set of *mutually exclusive* and *jointly comprehensive* subsets.
In this application, we are considering the population lineages (i.e., the tips of the population tree we give as input to |delineate|) as the elements of the initial set.
We want to organize this set of populations into distinct subsets, where each subset corresponds to a distinct species.

For example, imagine that that we have a system consisting of three populations, A, B, and C, which can be represented by the following set of three elements::

    {A, B, C}

What are the different "species delimitations" are possible with this system?
That is, what are the different ways we can we organize this set of populations in to species?

    1.  One way would be to lump all the populations to the same species, i.e.::

        {{A,B,C}}

    2.  Another way might be to split off each population to its own distinct species, i.e.::

        {{A},{B},{C}}

    3.  Another way might be to have A and B assigned to the same species, but C to a different one::

        {{A,B},{C}}

    4.  Or, alternatively, A and C to the same species, but B to a different one::

        {{A,C},{B}}

    5.  Or, finally, B and C to the same species, but A to a different one::

        {{A},{B,C}}

Each of these five arrangements is a distinct partition of the original set of three elements.
Thus, each possible partition of set of population lineages represents a different way of organizing that set of population lineages into species.
|delineate| will report the probabilities of each of these partitions for your data.
The partition with the highest probability represents the maximum likelihood estimate (MLE) of the species identities of the populations in your system.

For a given set of lineages, there are as many possible partitions as there are distinct groupings of lineages into subsets.
Each subset of lineages in a particular partition represents a distinct species, and, conversely, the membership of a particular lineage in a particular subset of a given partition is an assignment of that lineage to the species identity represented by that subset.
Thus, each partition represents a particular delimitation of lineages into species.
The boundaries between species then simply correspond to the boundaries between the subsets of the preferred partition.
The goal of our analysis is to identify the partition that best fits our data.

..
    The number of partitions possible for a set increases with the number of elements in the set.
    In fact, it increases very, very, very, very, very, very, `*very* rapidly <https://mathworld.wolfram.com/BellNumber.html>`_.

The Diversification Process
===========================

How are the probabilities of these different partitions of population lineages calculated?
To understand that we have to understand how |delineate| models the diversification process.
|delineate| models diversification with the following events:

-   *population splitting*
-   *population extinction*
-   *evolution of reproductive isolation of an independent population lineage*

Population splitting, i.e. the fragmentation of an ancestral population lineage into two independent daughter lineages, corresponds to a birth event on a tree of populations (from one or more species) .
Each of the two daughter lineages proceeds on its own independent trajectory, still part of the same species as the parent population lineages, until they either themselves split into further daughter lineages, go extinct (corresponding to death events on a tree of populations), or develop into a distinct species (from either parent or sister population).
As modeled by |delineate| the formation of new independent population lineages (birth events on a tree), the loss of these lineages due to extinction (death events on a tree), and the development of a lineage into a full independent species, are modeled as three independent processes, each with their own distinct rate.
Birth events on a tree can be considered the (potential) origin of a new species, and the the development of full independent species status as the completion of this phenomenon.
As such, with |delineate| speciation is not an instantaneous event (such as that modeled by a simple birth-death branching model underlying a species tree), but an extended process, starting with the lineeage splitting and ending with the speciation completion event.

Formation of New Species
========================

We can consider the speciation completion process to "play out" over a birth-death tree, much in the same way as character evolution might occur on a phylogeny.
Seen in this way, the tree of population lineages grows through a standard birth-death process, with a fixed birth rate and death rate representing the rate of population lineage splitting and extinction  respectively.
As noted above, a splitting event on this phylogeny does *not* represent a speciation event --- it just represents an ancestral population fragmenting into two isolated daughter populations.
Both daughter populations (as well as the nominally extinct ancestor population) belong to the same species.
The splitting event can be considered the initiation of potential speciation, however, as the isolation of the two daughter lineages essentially provides an opportunity for one or both of them to develop full reproductive isolation and thus achieve full independent species status.

Development of full independent species status is when a *speciation completion event* occurs on one of the daughter lineages before it goes extinct or itself split.
Seen in this way, a "species" in |delineate| is a set of population lineages on the population tree in which there is *no* speciation completion event on the edges connecting them on the tree.
If there is at least one speciation completion event on the edge path between two lineages, then the two lineage are in different species.

.. figure:: images/diversification.png
    :alt: The diversification process as modeled by |delineate|.

    Lineage splitting events correspond to the formation of new population lineages, not species, through restrictions in gene flow in an ancestral population (V1).
    These lineages may themselves give rise to other population lineages (V2 through V9), or go extinct (X1 through X3).
    Population lineages develop into an independent species at a fixed background rate, providing they are not otherwise lost  (i.e., there is duration between the initiation and completion of speciation).
    Changes in status from incipient to full or good species are marked by speciation completion events, shown by the blue bars.
    A "species" is thus made up of one or more population lineages not separated from one another by a speciation event.
    In this example, five speciation completion events divide the seven extant populations into four species: {A,B}, {C}, {D,E}, and {F,G}.

The Speciation Completion Rate
==============================

Speciation completion, i.e. the transition of an incipient species to full species status, completing a trajectory that started with its original splitting from a parent or sister population, proceeds at a rate given by the *speciation completion rate*.
This rate is one of the critical parameters that inform the probability of different species partitions, i.e. the different possible combinations of assignments of species identities to the various population lineages in the system.
For example, with a high species completion rate, partitions with more species would be more probable than partitions with fewer species as we would expect there to be more speciation completion events to have occured on the tree.
Conversely, with a low species completion rate, partitions with fewer species would be more probably than partitions with more species.

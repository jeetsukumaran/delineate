####################
Fundamental Concepts
####################

"Partitions"
============

The primary output of |delineate| is essentially calculation of the probabilities of different *partitions* of the population lineages in the study.
The term "partition" here comes from basic set theory: it is a division of a set of elements into a set of *mutually exclusive* and *jointly comprehensive* subsets.
In this application, we are considering the population lineages (i.e., the tips of the population tree we give as input to |delineate|) as the elements of the initial set.
We want to organize this set of populations into distinct subsets, where each subset corresponds to a distinct species.

For example, imagine that that we have a system consisting of three populations: A, B, and C.
How many different "species delimitations" are possible with this system?
That is, how many different ways can we organize this set of populations in to species?
One way would be to assign all the populations to the same species, i.e.::

{{A,B,C}}

Another way might be to assign each population to its own distinct species, i.e.::

{{A},{B},{C}}

Another way might be to have A and B assigned to the same species, but C to a different one::

{{A,B},{C}}

Or, alternatively, A and C to the same species, but B to a different one::

{{A,C},{B}}

Or, finally, B and C to the same species, but A to a different one::

{{A},{B,C}}

Each of these five arrangements is a distinct partition of the original set of three elements.
Thus, each possible partition of set of population lineages represents a different way of organizing that set of population lineages into species.
|delineate| will report the probabilities of each of these partitions for your data.
The partition with the highest probability represents the maximum likelihood estimate (MLE) of the species identities of the populations in your system.

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

We can imagine the speciation completion process to "play out" over a birth-death tree, much in the same way as character evolution might occur on a phylogeny.



Sampling Design
===============

|delineate| requires a fundamentally different way to thinking how we sample data for species delimitation studies.

|delineate| ideally should be provided with data that *includes samples from as many populations as possible across the system* being studied. This is in contrast from standard practice from other approaches, in which typically one or two examplar population sample per putative species are included. The theoretical ideal would be to include *every* population of all species in the system, known or unknown, i.e., to capture all population isolation or splitting events. Of course, we do not expect to achieve this theoretical ideal in practice, but it is certainly something to aspire to. The key point is restricting our population/species sampling to a few examplar populations per species is something we want to move away from.

|delineate| also requires that we have know the species identities of at least *some* of our population lineages. This is, again, in contrast to other approaches to species delimitation, which might be quite happy analyzing an entire data set with no known fixed species identities. With |delineate| we should design our sampling scheme to including a much broader range of species than just the ones we are interested in delimiting, and should include populations belong to species in which we are quite confident regarding their species identities. These other species --- or, to be more precise, population lineages for which the species identities are known --- are critical to allowing |delineate| to "learn" about the speciation process.

From Individuals to Populations to Species
==========================================

The |delineate| package itself represents the final step in an analytical pipeline.
Starting with a collection of genetic alignments representing multiple genes sampled from multiple individuals from multiple populations, each step of the pipeline groups the data into successively higher levels of organization, from populations to species.

A typical species delimitation analysis would consist of the following three steps:

    1.  **Identification (Delimitation) of *Population* Units:** First, we would carry out a |BPP|_ analysis to identify population units by aggregating individuals into populations under the multipopulation coalescent model. We would typically hope to sample at least a few genes from two to ten individuals per population, with multiple populations per putative species.
        We would then use |BPP|_ to organize these individuals into populations.
        Note that |BPP|_ terminology uses the term "species" and "populations" interchangeably. This can be confusing, but it is important to keep this in mind.

    2.  **Organization of the Population Units into a Population Phylogeny:** Then, we would carry out a |StarBeast2|_ analysis using the groupings identified by |BPP| as "species", to estimate an ultrametric phylogeny with those groupings as tips (i.e, a population phylogeny).
        Once we have decided what our population units are, we will use |StarBeast2|_ to infer an ultrametric population tree to use as input. Here, again, while |StarBeast2| uses the terminology "species" to reference to groupings of individuals, we should bear in mind that we are still dealing with population. We will use the units identified as populations by BPP as the "species" grouping in |StarBeast2|.

    3.  **Calculating the Probability of Species Assignments:** Finally, the actual species delimitation analysis itself: a |delineate| analysis to calculate the probabilities of different groupings of population tips of the population tree into species.
        The population tree resulting from |StarBeast2|_ forms the one of the mandatory inputs for |delineate|. The species identities for the subset of population lineages for which these are known forms the other. Running |delineate| will then report on the probabilities of different species assignments for the remaining lineages, i.e. for the ones for which we do not know or specifies the species identities.


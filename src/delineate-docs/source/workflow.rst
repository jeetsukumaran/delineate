########
Workflow
########

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



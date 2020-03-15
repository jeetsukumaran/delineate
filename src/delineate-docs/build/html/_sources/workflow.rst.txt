#########################
Study Design and Workflow
#########################

1. Sampling Design

   DELINEATE requires a fundamentally different way to thinking how we
   sample data for species delimitation studies.

   DELINEATE ideally should be provided with data that *includes samples
   from as many populations as possible across the system* being
   studied. This is in contrast from standard practice from other
   approaches, in which typically one or two examplar population sample
   per putative species are included. The theoretical ideal would be to
   include *every* population of all species in the system, known or
   unknown, i.e., to capture all population isolation or splitting
   events. Of course, we do not expect to achieve this theoretical ideal
   in practice, but it is certainly something to aspire to. The key
   point is restricting our population/species sampling to a few
   examplar populations per species is something we want to move away
   from.

   DELINEATE also requires that we have know the species identities of
   at least *some* of our population lineages. This is, again, in
   contrast to other approaches to species delimitation, which might be
   quite happy analyzing an entire data set with no known fixed species
   identities. With DELINEATE we should design our sampling scheme to
   including a much broader range of species than just the ones we are
   interested in delimiting, and should include populations belong to
   species in which we are quite confident regarding their species
   identities. These other species --- or, to be more precise,
   population lineages for which the species identities are known ---
   are critical to allowing DELINEATE to "learn" about the speciation
   process.

2. From Individuals to a Populations

   We would typically hope to sample at least a few genes from two to
   ten individuals per population per species. We would then use
   `BPP <https://github.com/bpp/bpp>`__ to organize these individuals
   into populations. Note that `BPP <https://github.com/bpp/bpp>`__
   terminology uses the term "species" and "populations"
   interchangeably. This can be confusing, but it is important to keep
   this in mind.

3. Estimating the Population Tree

   Once we have decided what our population units are, we will use
   `StarBeast2 <https://taming-the-beast.org/tutorials/starbeast2-tutorial/>`__
   to infer an ultrametric population tree to use as input. Here, again,
   while
   `StarBeast2 <https://taming-the-beast.org/tutorials/starbeast2-tutorial/>`__
   uses the terminology "species" to reference to groupings of
   individuals, we should bear in mind that we are still dealing with
   population. We will use the units identified as populations by BPP as
   the "species" grouping in [StartBeast2].

4. Estimating the Probability of Species Assignments

   The population tree resulting from
   `StarBeast2 <https://taming-the-beast.org/tutorials/starbeast2-tutorial/>`__
   forms the one of the mandatory inputs for DELINEATE. The species
   identities for the subset of population lineages for which these are
   known forms the other. Running DELINEATE will then report on the
   probabilities of different species assignments for the remaining
   lineages, i.e. for the ones for which we do not know or specifies the
   species identities.

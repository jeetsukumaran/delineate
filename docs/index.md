# DELINEATE: Species Delimitation

## Introduction

DELINEATE is a program that implements the DELINEATE approach to species delimitation (Sukumaran, Holder, and Knowles, 2019).
This approach integrates an explicit model of speciation into the "censored" or "multispecies" coalescent model to organize a set of population lineages sampled from one or more species into mutually-exclusive and jointly-comprehensive subsets, where each subset of population lineages represents a distinct species.

## Requirements

-   [Python] (3.7 or greater)
-   [DendroPy] (4.4.0 or greater)

    If you already have [Python] 3 available, you can install or upgrade [DendroPy] by:

~~~
python3 -m pip install -U git+https://github.com/jeetsukumaran/DendroPy.git
~~~

## Installation or Upgrade

You can run the following command to install or upgrade DELINEATE directly from the package repositories::

    python3 -m pip install -U git+https://github.com/jeetsukumaran/delineate.git

## Data Requirements

DELINEATE requires two items of data:

-   A *population tree*
-   A *species assignment table*

### Input Data: Population Tree

This is rooted ultrametric tree where each tip lineage represents a population or deme.
This tree is typically obtained through a classical "censored" or "multispecies" coalescent analysis, such as results from BP&P (either mode A01 or A10) or StarBeast.
The tree can be be specified either in NEXUS or Newick format.

### Input Data: Species Assignment Table

This is a tab-delimited plain text file with at least three columns:

-   "lineage"
-   "species"
-   "status"

There may be more than these three columns, but these columns are mandatory (and all other columns will be ignored).
The order of columns does not matter, as the DELINEATE programs will use the labels specified in the header row (see below) to identify the columns.

The first row is the header row: i.e., column labels.
Subsequent rows will map each tip in the population tree (the "lineage" column) to a species label ("species") as well as an indication of whether this specis assignment is known or not ("status").
*Every* tip in the population tree must be represented by a row (and no more than one row) in this table.
If species assignments are not known for some population lineages (as would be expected in species discovery type analyses), then the "species" field can be left blank or populated with an arbitrary value, such as a "?" or "NewSp.?" etc., but the "status" field should be set to "0".
Species assignments that *are* known will have the status field set to "1".

The following example illustrates the structure and semantics of the species assignment table:

| lineage | species | status |
|---------|---------|--------|
| able1   | R_able  | 1      |
| able2   | R_able  | 1      |
| able3   | R_able  | 1      |
| new1    | ???     | 0      |
| baker1  | R_baker | 1      |
| baker2  | R_baker | 1      |
| easy1   | R_easy? | 0      |
| easy2   | R_easy? | 0      |
| oboe1   | unknown | 0      |
| oboe2   | unknown | 0      |

In this case, the population tree should consist of exactly 9 tips with the following labels: "able1", "able2", "able3", "new1", "baker1", "baker2", "easy1", "easy2", "oboe1", "oboe2".
The "1's" in the "status" column indicate species assignments that should be taken as known *a priori* by the program.
The "0's" in the "status" column indicate populations of unknown species uncertainties: DELINEATE will assign species identities to these lineages (either existing ones, such as "R_able" or "R_baker", or establish entirely new species identities for them).
Note that the species labels are ignored for population lineages with a "0" status --- these are just there for user book-keeping or reference.

## Running a Species Delimitation Analysis

### Basic Run

Given a population lineage tree file, "population-tree.nex", and a species assignment table file "species-mappings.tsv", then the following command will run a DELINEATE analysis on the data::

~~~
delineate-estimate partitions --tree-file population-tree.nex --config-file data1.tsv
~~~

or, using the short-form options:

~~~
delineate-estimate partitions -t population-tree.nex -c data1.tsv
~~~

This command has the following components:

-   ``delineate-estimate``:
    This is the name of the program to be run.
-   ``estimate``:
    This is the command or operation that the program will be running.
-   ``--tree-file population-tree.nex`` or ``-t population-tree.nex``:
    The ``--tree-file`` flag, or its short-form synonym, ``-t``, specifies that the next element will be the path to tree file with data on the population lineage tree.
    In this example, the file is located in the current working directory, i.e., *population-tree.nex*.
    If it was in another directory, then the path could be */home/bilbo/projects/orc-species-delimitation/data1/population-tree.nex*, for example.
-   ``--config-file data1.tsv`` or ``-c data1.tsv``
    The ``--config-file`` flag, or its short-form synonym, ``-c``, specifies that the next element will be the path to species assignment configuration file.
    In this example, the file is located in the current working directory, i.e., *delineate-species.tsv*.
    Again, if it was in another directory, then the path could be */home/bilbo/projects/orc-species-delimitation/data1/data1.tsv*, for example.

### Basic Run Output

Executing this command will run the DELINEATE analysis and will produce the following output files:

-   "*data1.delimitation-results.json*"
-   "*data1.delimitation-results.trees*"

More generally, unless the ``-o`` or ``--output-prefix`` flag (see below) is used to explicity specify an alternate output prefix for all results generated, the files will take on a prefix given by the file stemname of the configuration file, "*data1*" in the this example.

The first file, with the general name of "*&lt;output-prefix&gt;.delimitation-results.json*", is the primary results file.
As can be inferred from its extension, it is a [JSON] format text file, and it consists of a single a dictionary.
The dictionary provides information on the estimated speciation completion rate as well as the probabilities of all the possible partitions of the population lineage leafset into species sets, given the species assignment constraints, ranked by the probability of each partition.

The second file, with the general name of "*&lt;output-prefix&gt;.delimitation-results.trees*",  provides supporting results.
Basically this is a collection of trees, with one tree for each partition considered.
The topology of the trees are identical, corresponding to the topology of the input tree (i.e., the population lineage tree), as are the tip labels.
However, the tips have associated with them some extra metadata that will be available for viewing in a program like [FigTree].
Most important of this is ``species``, i.e., the label corresponding to the identity of the species assignment in that partition.
In the case of species assignments that are constrained (i.e., status indicated by "1"), these will be identical to the assignment and invariant across all partitions, of course.
However, in the case of population lineages of *unknown* species affinities (i.e., status indicated by "0"), this may be an existing species label (if the population lineage was assigned to an existing species in the partition under consideration) or a new, arbitrary species label (if the population lineage was assigned to a new distinct species in the partition under consideration).
In addition, in [FigTree] you can also choose to have the branches colored by "status", and this will highlight population lineages of (*a priori*) known vs unknown species affinities, and thus quickly identify the assigned species identities of the lineages of interest.

### Calculating the Marginal Probability of Conspecificity for a Subset of Taxa

An analysis estimating the probabilities of different partition gives the the \textit{joint} probability for different organizations of population lineages into subsets, each constituting a distinct species.
We might be interested, instead, in the marginal probability of conspecificity of a few population lineages.
That is, we are asking the question, "What is the probability that the following population lineages are conspecific?"
To answer the question we simply have to sum up the probabilities of individual partitions in which those population lineages are found to be conspecific.
We provide an application to this, ``delineate-summarize``.
To run it, you simply need to provide it the JSON results file that is produced by a ``delineate-estimate`` analysis, followed by the list of taxa for which you want to calculate the marginal probability of conspecificity.
Note that if your taxon labels have spaces or special characters in them (tsk tsk), you need to wrap your labels in quotes.
E.g., assuming you have run a DELINEATE analysis using ``delineate-estimate``, and the analysis produced the two following files:

-   "*dyna1.delimitation-results.json*"
-   "*dyna1.delimitation-results.trees*"

Then, to calculate the marginal probability that, for e.g., the population lineages "DGRP1" and "DGRR1" are conspecific, you would run the following command:

~~~
delineate-summarize conspecificity -r dyna1.delimitation-results.json DGRP1 DGRR1
~~~

Or the marginal probability that the population lineages "DhtT9", "Dhy3Br", "Dhy6", and "Dhym5" are conspecific:

~~~
delineate-summarize conspecificity -r dyna1.delimitation-results.json DhtT9 Dhy3Br Dhy6 Dhym5
~~~

[Python]: https://www.python.org/
[DendroPy]: https://dendropy.org/
[JSON]: https://en.wikipedia.org/wiki/JSON
[FigTree]: http://tree.bio.ed.ac.uk/software/figtree/


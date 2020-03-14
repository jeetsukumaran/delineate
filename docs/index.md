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
You can also have the tip labels to displayed both the assigned species as well as the population label together ("species-lineage" or "lineage-species" options).
In the case of species assignments that are constrained (i.e., status indicated by "1"), these will be identical to the assignment and invariant across all partitions, of course.
However, in the case of population lineages of *unknown* species affinities (i.e., status indicated by "0"), this may be an existing species label (if the population lineage was assigned to an existing species in the partition under consideration) or a new, arbitrary species label (if the population lineage was assigned to a new distinct species in the partition under consideration).
In addition, in [FigTree] you can also choose to have the branches colored by "status", and this will highlight population lineages of (*a priori*) known vs unknown species affinities, and thus quickly identify the assigned species identities of the lineages of interest.

### More Options

All subcommands for ``delineate-estimate`` will be shown with ``--help`` option:

~~~
delineate-estimate --help
~~~

While special options for the ``partitions`` will be shown by:

~~~
delineate-estimate partitions --help
~~~

which results in:

~~~
usage: delineate-estimate partitions [-h] -t TREE_FILE [-c CONFIG_FILE]
                                     [-f {nexus,newick}]
                                     [--underscores-to-spaces] [-u]
                                     [--speciation-completion-rate-estimation-min #.##]
                                     [--speciation-completion-rate-estimation-max #.##]
                                     [--speciation-completion-rate-estimation-initial #.##]
                                     [-o OUTPUT_PREFIX]
                                     [--extra-info EXTRA_INFO_FIELD_VALUE]
                                     [-s SPECIATION_COMPLETION_RATE] [-P #.##]
                                     [--report-mle-only] [-p #.##] [-I]
                                     [--no-translate-tree-tokens]
                                     [-l {lineage,species,lineage-species,species-lineage,clear}]
                                     [--store-relabeled-trees {lineage-species,species-lineage}]

Given a known population tree and optionally a speciation completion rate,
calculate the probability of different partitions of population lineages into
species, with the partition of the highest probability corresponding to the
maximum likelihood species delimitation estimate.

Command Options:
  -h, --help            show this help message and exit
  --underscores-to-spaces, --no-preserve-underscores
                        Convert underscores to spaces in tree lineage labels (if
                        not proected by quotes). in the configuration file will
                        not be modified either way. You should ensure consistency
                        in labels between the tree file and the configuration
                        file.

Source Options:
  -t TREE_FILE, --tree-file TREE_FILE
                        Path to tree file.
  -c CONFIG_FILE, --config-file CONFIG_FILE
                        Path to configuration file.
  -f {nexus,newick}, --tree-format {nexus,newick}
                        Tree file data format (default='nexus').

Estimation Options:
  -u, --underflow-protection
                        Try to protect against underflow by using special number
                        handling classes (slow).
  --speciation-completion-rate-estimation-min #.##, --smin #.##
                        If estimating speciation completion rate, minimum
                        boundary for optimization window [default: 1e-08].
  --speciation-completion-rate-estimation-max #.##, --smax #.##
                        If estimating speciation completion rate, maximum
                        boundary for optimization window [default: 10 x lineage
                        tree pure birth rate].
  --speciation-completion-rate-estimation-initial #.##, --sinit #.##
                        If estimating speciation completion rate, initial value
                        for optimizer [default: 0.01 x lineage tree pure birth
                        rate].

Output Options:
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix for output file(s).
  --extra-info EXTRA_INFO_FIELD_VALUE
                        Extra information to append to output (in format <FIELD-
                        NAME>:value;)
  -I, --tree-info       Output additional information about the tree.

Model Options:
  -s SPECIATION_COMPLETION_RATE, --speciation-completion-rate SPECIATION_COMPLETION_RATE
                        The speciation completion rate. If specified, then the
                        speciation completion rate will *not* be estimated, but
                        fixed to this.

Report Options:
  -P #.##, --report-cumulative-probability-threshold #.##
                        Do not report on partitions outside of this constrained
                        (conditional) cumulative probability.
  --report-mle-only     Only report maximum likelihood estimate.
  -p #.##, --report-probability-threshold #.##
                        Do not report on partitions with individual constrained
                        (conditional) probability below this threshold.

Tree Output Options:
  --no-translate-tree-tokens
                        Write tree statements using full taxon names rather than
                        numerical indices.
  -l {lineage,species,lineage-species,species-lineage,clear},
  --figtree-display-label {lineage,species,lineage-species,species-lineage,clear}
                        Default label to display in FigTree for the trees in the
                        main result file.
  --store-relabeled-trees {lineage-species,species-lineage}
                        Create an additional result file, where trees are written
                        with their labels actually changed to the option selected
                        here.
~~~

Some of these options are explained in detail below:

*   ``-o`` or ``--output-prefix``

    By default ``delineate-estimate`` will create output files in the current
    (working) directory with a filename stem given by the filename stem of the
    configuration file. If you out want to specify a different path and/or
    directory, you would use this option. For e.g.,

    ~~~
    delineate-estimate partitions \
        -t poptree.nex \
        -c data1.tsv \
        -o /arrakis/workspace/results/run1
    ~~~

    or

    ~~~
    delineate-estimate partitions \
        --tree-file poptree.nex \
        --config-file data1.tsv \
        --output-preifx /arrakis/workspace/results/run1
    ~~~

    will result in output being written to the following locations:

    -   "*/arrakis/workspace/results/run1.delimitation-results.json*"
    -   "*/arrakis/workspace/results/run1.delimitation-results.trees*"

*   ``-P`` or ``report-probability-threshold``

    By default, ``delineate-estimate`` will enumerate *all* possible partitions and their probabilities.
    For many studies, this will result in massive datafiles, hundreds to thousands of gigabytes, as millions or millions of millions or more partitions are written out.
    Remember, the number of partitions for a dataset of $n$ lineages is the $n^{th}$ Bell number.
    For most studies, the vast majority of these partitions will be of very, very, very low probability.
    Considerable savings in time, disk space, and a significant slow down of the heat death of the universe can be achieved by restricting the report to only most probable partitions that collectively contribute to most of the probability.
    This option allows you to effect this restriction and achieve these saving.
    To restrict the report to the most probable partitions that collectively result in a cumulative probability of 0.95, for e.g.,

    ~~~
    delineate-estimate partitions -t poptree.nex -c data1.tsv -P 0.95
    ~~~

    or

    ~~~
    delineate-estimate partitions --tree-file poptree.nex --config-file data1.tsv --report-cumulative-probability-threshold 0.95
    ~~~

    Note that one issue that might result from this restriction is that, when summing up the marginal probabilities of the statues of a particular lineage or collection of lineages (new species, conspecificity, etc.), exclusion of these lower probability partitions may have an effect on the accuracy of the marginal probabilties.
    However, in most cases this is probably not going to be very consequential, especially with a sufficiently high threshold such that there is very little probability mass remaining in the unreported part of the parameter space.

*   ``-u`` or ``--underflow-protection``

    For some data sets, especially largers ones with relatively few constraints, ``delineate-estimate`` may crash due to underflow errors: the probabilities of some of the partitions are simply too low to be handled correctly by our silicon counterparts (and log-transforming is not possible as we are summing some of this sub-probabilities).
    The wonderful [Python] Standard Library provides a really elegant solution to handling all sorts of numerical woes including this -- a specialized number class ([Decimal]).
    While I have not done extensive benchmarking, however, I am pretty sure that usage of this class instead of the ``float`` primitive brings with it a performance cost.
    As such, ``delineate-estimate`` does not use [Decimal] by default.
    However, if you *do* find your analysis crashing due to underflow errors, then you may have to take the performance hit to see the journey's end.
    In these cases, specifying the underflow protection flag will (hopefully) get you to that end through a more robust albeit slower road.

    ~~~
    delineate-estimate partitions -t poptree.nex -c data1.tsv -u
    ~~~

    or

    ~~~
    delineate-estimate partitions --tree-file poptree.nex --config-file data1.tsv --underflow-protection
    ~~~

## Summarizing the Marginal Probability of Conspecificity and New Species Status for a Subset of Taxa

An analysis estimating the probabilities of different partition gives the \textit{joint} probability for different organizations of population lineages into subsets, with each subset constituting a distinct species.
The maximum likelihood estimate of the species delimitation is given by the partition with the highest probability, and this is easily available from the JSON results file or the tree file (it is the first partition in the JSON file, and the first tree in the tree file).

However, we might be interested instead in a number of other questions, such as:

*   What the marginal probability of conspecificity of a particular set of population lineages? That is, we are asking the question, "What is the probability that the following population lineages are conspecific?"
*   Or we might be interested in a stricter condition, i.e. *exclusive conspecificity*, i.e., the marginal probability that a group of lineages form an exclusive complete species unto themselves, including no other lineages.
*   Or we might like to know the marginal probability that one or more lineages might be assigned to new species, not provided by us as known identity in the constraints to the analysis.
*   Or we might like to know the marginal probability that one or more lineages might form a distinct, self-contained, *exclusive* AND *new* species (as opposed to just an exclusive species, or just a non-exclusive part of new species)

These questions are readily answered by summing the probabilities of various partitions that meet the above conditions as given in the results file.
We provide an application to this, ``delineate-summarize``.
To run it, you simply need to provide it the JSON results file that is produced by a ``delineate-estimate`` analysis, followed by the list of lineages for which you want to calculate the marginal probabilities for one or more of the above conditions.
Note that if your taxon labels have spaces or special characters in them (tsk tsk), you need to wrap your labels in quotes.
E.g., assuming you have run a DELINEATE analysis using ``delineate-estimate``, and the analysis produced the two following files:

-   "*biologicalconcept.delimitation-results.json*"
-   "*biologicalconcept.delimitation-results.trees*"

Then, to calculate the marginal probability that, for e.g., the population lineages "DhrM1", "DHHG2", and "DHHG2" are conspecific, exclusively conspecfic, are part of a new species, or form an exclusive new species you would run the following command:

~~~
delineate-summarize -r biologicalconcept.conf.delimitation-results.json DhrM1 DHHG2 DhhD1
~~~

which might result in something like this:
~~~
[delineate-summarize] 25 lineages defined in results file: 'CCVO1', 'DGRP1', 'DGRR1', 'DHHG2', 'DhbV2', 'DheCO6', 'DheP8', 'DhhD1', 'DhlCO1', 'DhlP7', 'DhmB2Br', 'DhpB1Br', 'DhrM1', 'DhrSL5', 'DhsG1', 'DhsH3', 'DhsP1', 'DhtT9', 'Dhy3Br', 'Dhy6', 'Dhym5', 'Dma2', 'DtTN1', 'ERVL2', 'YSNE1'
[delineate-summarize] 5 species defined in configuration constraints, with 12 lineages assigned:
    [  1/5  ] 'ecu'    (3 lineages)
    [  2/5  ] 'hy'     (2 lineages)
    [  3/5  ] 'lic'    (3 lineages)
    [  4/5  ] 'ma'     (1 lineages)
    [  5/5  ] 'sep'    (3 lineages)
[delineate-summarize] 12 out of 25 lineages assigned by constraints to 5 species:
    [  1/12 ] 'DheCO6'    (SPECIES: 'ecu')
    [  2/12 ] 'DheP8'     (SPECIES: 'ecu')
    [  3/12 ] 'YSNE1'     (SPECIES: 'ecu')
    [  4/12 ] 'Dhy3Br'    (SPECIES: 'hy')
    [  5/12 ] 'Dhy6'      (SPECIES: 'hy')
    [  6/12 ] 'DhlCO1'    (SPECIES: 'lic')
    [  7/12 ] 'DhlP7'     (SPECIES: 'lic')
    [  8/12 ] 'ERVL2'     (SPECIES: 'lic')
    [  9/12 ] 'Dma2'      (SPECIES: 'ma')
    [ 10/12 ] 'DhsG1'     (SPECIES: 'sep')
    [ 11/12 ] 'DhsH3'     (SPECIES: 'sep')
    [ 12/12 ] 'DhsP1'     (SPECIES: 'sep')
[delineate-summarize] 13 out of 25 lineages not constrained by species assignments:
    [  1/13 ] 'CCVO1'
    [  2/13 ] 'DGRP1'
    [  3/13 ] 'DGRR1'
    [  4/13 ] 'DHHG2'
    [  5/13 ] 'DhbV2'
    [  6/13 ] 'DhhD1'
    [  7/13 ] 'DhmB2Br'
    [  8/13 ] 'DhpB1Br'
    [  9/13 ] 'DhrM1'
    [ 10/13 ] 'DhrSL5'
    [ 11/13 ] 'DhtT9'
    [ 12/13 ] 'Dhym5'
    [ 13/13 ] 'DtTN1'
[delineate-summarize] 80271 partitions found in results file, with total constrained cumulative probability of 0.9500004131628125
[delineate-summarize] Reading focal lineages from arguments
[delineate-summarize] 3 focal lineages defined: 'DHHG2', 'DhhD1', 'DhrM1'
[delineate-summarize] 52642 out of 80271 partitions found with focal lineages conspecific
[delineate-summarize] Marginal constrained probability of focal lineages conspecificity: 0.8314956644192
[delineate-summarize] Marginal constrained probability of focal lineages *exclusive* conspecificity: 0.05232352098265284
[delineate-summarize] Marginal constrained probability of focal lineages being collectively an exclusive new species: 0.05232352098265284
[delineate-summarize] Marginal constrained probability of focal lineages being collectively *part* (i.e., non-exclusively) of a new species: 0.3777092230039129
[delineate-summarize] Marginal constrained probability of focal lineages being collectively *part* (i.e., non-exclusively) of a predefined species: 0.4537864414152871
[delineate-summarize] WARNING: cumulative constrained probability in results file is only 0.9500004131628125. Not all partitions might have been included, and probability summarizations reported should not be considered as accurate.
{
    "lineages": ["DHHG2", "DhhD1", "DhrM1"],
    "marginal_probability_of_conspecificity": 0.8314956644192,
    "marginal_probability_of_exclusive_conspecificity": 0.05232352098265284,
    "marginal_probability_of_new_species": 0.3777092230039129,
    "marginal_probability_of_existing_species": 0.4537864414152871,
    "marginal_probability_of_exclusive_new_species": 0.05232352098265284
}
~~~

As can be seen, the complete report includes details on the configuration constraints etc. as well as the various summarized marginal probabilities.
The final result is written (by default) in JSON format to the standard output.
The output format and details can be changed by specifying different options to the ``delineate-summarize`` program. For a full list of these, type:

~~~
delineate-summarize --help
~~~

Note that the results also report the marginal probabilities that the set of taxa constitute collectively part of a new species (i.e., a species definition not provided to DELINEATE as part of the constraints). So we can also just pass in the name of a single population lineage to see the marginal probability it was placed in a species distinct from any named in the contraints:

~~~
delineate-summarize -r biologicalconcept.delimitation-results.json DhtT9
~~~

[Python]: https://www.python.org/
[DendroPy]: https://dendropy.org/
[JSON]: https://en.wikipedia.org/wiki/JSON
[FigTree]: http://tree.bio.ed.ac.uk/software/figtree/
[Decimal]: https://docs.python.org/3.8/library/decimal.html#decimal-objects

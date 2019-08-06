DELINEATE: Species Delimitation
===============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Introduction
------------

DELINEATE is a program that implements the DELINEATE approach to species delimitation (Sukumaran, Holder, and Knowles, 2019).
This approach integrates an explicit model of speciation into the "censored" or "multispecies" coalescent model to organize a set of population lineages sampled from one or more species into mutually-exclusive and jointly-comprehensive subsets, where each subset of population lineages represents a distinct species.

Requirements
------------

*   Python_ (3.7 or greater)
*   DendroPy_ (4.4.0 or greater)

    If you already have Python 3 available, you can install DendroPy_ by::

        python3 -m pip install git+https://github.com/jeetsukumaran/DendroPy.git

Installation
------------

You can run the following command to install DELINEATE directly from the package repositories::

    python3 -m pip install git+https://github.com/jeetsukumaran/delineate.git

Input Data
----------

DELINEATE requires two items of data:

*   A *population tree*
*   A *species assignment table*

Input Data: Population Tree
...........................

This is rooted ultrametric tree where each tip lineage represents a population or deme.
This tree is typically obtained through a classical "censored" or "multispecies" coalescent analysis, such as results from BP&P (either mode A01 or A10) or StarBeast.
The tree can be be specified either in NEXUS or Newick format.

Input Data: Species Assignment Table
....................................

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

+---------+---------+--------+
| lineage | species | status |
+=========+=========+========+
| able1   | R_able  | 1      |
+---------+---------+--------+
| able2   | R_able  | 1      |
+---------+---------+--------+
| able3   | R_able  | 1      |
+---------+---------+--------+
| new1    | ???     | 0      |
+---------+---------+--------+
| baker1  | R_baker | 1      |
+---------+---------+--------+
| baker2  | R_baker | 1      |
+---------+---------+--------+
| easy1   | R_easy? | 0      |
+---------+---------+--------+
| easy2   | R_easy? | 0      |
+---------+---------+--------+
| oboe1   | unknown | 0      |
+---------+---------+--------+
| oboe2   | unknown | 0      |
+---------+---------+--------+

In this case, the population tree should consist of exactly 9 tips with the following labels: "able1", "able2", "able3", "new1", "baker1", "baker2", "easy1", "easy2", "oboe1", "oboe2".
The "1's" in the "status" column indicate species assignments that should be taken as known *a priori* by the program.
The "0's" in the "status" column indicate populations of unknown species uncertainties: DELINEATE will assign species identities to these lineages (either existing ones, such as "R_able" or "R_baker", or establish entirely new species identities for them).

Running the Program
-------------------

Output
------

.. _Python: https://www.python.org/
.. _DendroPy: https://dendropy.org/

.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

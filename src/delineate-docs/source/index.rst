%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DELINEATE: Species Delimitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

|delineate| is a package that implements the |delineate| approach to species delimitation (Sukumaran, Holder, and Knowles, 2020). This approach integrates an explicit model of speciation into the multipopulation coalescent (also known as the "censored" or "multispecies" coalescent) model, to organize a set of population lineages sampled from one or more species into mutually-exclusive and jointly-comprehensive subsets, where each subset of population lineages represents a distinct species.

|delineate| requires |Python|_ 3.7 or greater to run.
Package and dependency management is greatly simplified if |pip|_ is installed as well.

Installation
============

Providing you have all the system requirements (see :doc:`the installation guide </install>`), you can use the following command to install or upgrade |delineate| directly from the package repositories::

    python3 -m pip install -U git+https://github.com/jeetsukumaran/delineate.git

Documentation
=============


.. toctree::
    :maxdepth: 3

    install.rst
    concepts.rst
    quickstart.rst
    workflow1.rst




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

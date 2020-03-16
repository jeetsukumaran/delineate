############
Installation
############

Installation or Upgrade
=======================

Providing you have all the system requirements (see below), you can use the following command to install or upgrade |delineate| directly from the package repositories::

    python3 -m pip install -U git+https://github.com/jeetsukumaran/delineate.git

System Requirements
===================

-  `Python <https://www.python.org/>`__ (3.7 or greater)

   Make sure you have `Python 3.7 or
   higher <https://www.python.org/downloads/>`__ available on your
   system. You can check this by typing the following in your operating
   system shell:

   ::

       $ python3

   If you get a response indicating that the operating system could not
   find or understand that command, or if it shows that the version is
   Python 3.6 or lesser, then please install an appropriate version of
   [Python3] from: https://www.python.org/downloads/ .

   (**NOTE**: I highly recommend that you consider using
   |Anaconda|_ instead of the system
   |Python|_ installation. As a general best
   practice, you should keep your operating system Python environments
   isolated from yours. Admittedly, this does require some extra effort
   and learning on your part (not just how to use Python, but how to
   manage Python environments using the
   `Anaconda <https://www.anaconda.com>`__ system). But it *really* is
   worth it! Please visit the `Anaconda <https://www.anaconda.com>`__
   site for details on how to setup and use
   `Anaconda <https://www.anaconda.com>`__ on your system.)

-  `pip <https://docs.python.org/3/installing/index.html>`__ (for Python
   3, i.e. ``pip3``)

   `pip <https://docs.python.org/3/installing/index.html>`__ should be
   installed by default along with your
   `Python <https://www.python.org/>`__ for all newer Python versions (>
   3.4), except, frustratingly enough, in the case of operating system
   Python's for some flavors of Linux. As I note above, I really
   encourage *all* users to use `Anaconda <https://www.anaconda.com>`__
   for their Python usage rather than the operating system Python. But
   if you *do* wish to use the operating system Python, then you will
   have to use your operating system package manager to install
   `pip <https://docs.python.org/3/installing/index.html>`__:

   ::

       sudo apt install python3-pip

   Note that, rather confusingly, Python 2 has its own version of
   `pip <https://docs.python.org/3/installing/index.html>`__, so if your
   operating system has Python 2 in addition to Python 3 side-by-side,
   you should ensure that you specifically have and will be running the
   Python 3 version by using the full version-qualified command name of
   ``pip3`` rather than just ``pip``. Type ``pip3 --version`` in your
   operation system shell to check:

   ::

       pip3 --version

Library Dependencies
====================

|delineate| makes use of the following Python libraries.

-  `DendroPy <https://dendropy.org/>`__ (4.4.0 or greater)
-  `SciPy <https://www.scipy.org/>`__ (1.18.1 or greater)
-  `NumPy <https://numpy.org/>`__ (1.4.1 or greater)

These dependencies should be installed automatically by running the
``pip`` command below. In case that does not work, you can manually
install these dependencies yourself by running:

::

    python3 -m pip install -U \
        scipy \
        numpy \
        git+https://github.com/jeetsukumaran/DendroPy.git

#! /usr/bin/env python

#from distutils.core import setup
from setuptools import setup
from delineate import __version__, __project__

setup(
    name=__project__,
    version=__version__,
    author="Jeet Sukumaran and Mark T. Holder",
    author_email="jeetsukumaran@gmail.com and mtholder@ku.edu",
    packages=["delineate", "delineate.test"],
    # scripts=["bin/delineate.py",],
    test_suite = "delineate.test",
    url="http://pypi.python.org/pypi/delineate/",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.txt").read(),
    # install_requires=[ ],
)

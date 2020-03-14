#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re
from setuptools import setup, find_packages

def _read(path_components, **kwargs):
    path = os.path.join(os.path.dirname(__file__), *path_components)
    if sys.version_info.major < 3:
        return open(path, "rU").read()
    else:
        with open(path, encoding=kwargs.get("encoding", "utf8")) as src:
            s = src.read()
        return s

project_init = _read(["src", "delineate", "__init__.py"])
__version__ = re.match(r".*^__version__\s*=\s*['\"](.*?)['\"]\s*$.*", project_init, re.S | re.M).group(1)
__project__ = re.match(r".*^__project__\s*=\s*['\"](.*?)['\"]\s*$.*", project_init, re.S | re.M).group(1)

setup(
    name=__project__,
    version=__version__,
    author="Jeet Sukumaran and Mark T. Holder",
    author_email="jeetsukumaran@gmail.com and mtholder@ku.edu",
    packages=find_packages("src"),
    package_dir={"": "src"},
    scripts=[
        "bin/delineate-check",
        "bin/delineate-estimate",
        "bin/delineate-summarize",
        ],
    test_suite = "tests",
    url="https://github.com/jeetsukumaran/delineate/",
    license="LICENSE.txt",
    description="Model-based species delimitation.",
    long_description=_read(["README.txt"]),
    install_requires=[
        "numpy>=1.18.1",
        "scipy>=1.4.1",
        #"DendroPy @ git+https://github.com:jeetsukumaran/DendroPy.git#egg=DendroPy",
        "DendroPy @ git+https://github.com:jeetsukumaran/DendroPy.git@development-master#egg=DendroPy",
        ]
)

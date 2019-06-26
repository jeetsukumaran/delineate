#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import argparse
from delineate import utility
from delineate import model

"""
Generate normalized JSON/delimited format configuration file from a JSON/delimited data file.
"""

__prog__ = os.path.basename(__file__)
__version__ = "1.0.0"
__description__ = __doc__
__author__ = 'Jeet Sukumaran and Mark T. Holder'
__copyright__ = 'Copyright (C) 2019 Jeet Sukumaran and Mark T. Holder.'

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("-c", "--config-file",
            action="store",
            help="Path to configuration file; can be JSON for columnar tab/comma-delimited file.")
    parser.add_argument("-d", "--config-file-delimiter",
            action="store",
            dest="delimiter",
            default=None,
            help="[Assuming non-JSON columnar file] configuration data delimiting character.")
    parser.add_argument("-t", "--tree_file",
            default=None,
            help="Path to tree file (optional; for validation)")
    parser.add_argument("-f", "--tree-format",
            dest="tree_format",
            type=str,
            default="nexus",
            choices=["nexus", "newick"],
            help="Tree file data format (default='%(default)s').")
    parser.add_argument("--output-format",
            choices=["json", "tsv", "csv"],
            default="json",
            help="Format for output file (default='%(default)s').")
    parser.add_argument("-o", "--output-prefix",
            action="store",
            default=None,
            help="Prefix for output file.")
    parser.add_argument("--ignore-extra-tree-lineages",
            action="store_true",
            default=False,
            help="Do not complain if not all tree lineages are specified in the configuration file.")
    parser.add_argument("--ignore-extra-configuration-lineages",
            action="store_true",
            default=False,
            help="Do not complain if not all configuration lineages are found on the tree.")
    args = parser.parse_args()
    if args.delimiter is not None and (
            args.delimiter.upper() == "TAB"
            or args.delimiter == '\\t'
            or args.delimiter == '\t'
            ):
        args.delimiter = "\t"
    logger = utility.RunLogger(name="delineate-configure",
            is_include_name=True,
            is_include_timestamp=False,
            log_to_stderr=True,
            stderr_logging_level=utility.logging.INFO,
            log_to_file=False,
            file_logging_level=utility.logging.INFO,
            )
    config_d = utility.parse_configuration(
            args=args,
            delimiter=args.delimiter,
            logger=logger)
    if args.output_format == "json":
        suffix = ".json"
        output_format = "json"
        output_delimiter="\t"
    elif args.output_format == "tsv":
        suffix = ".tsv"
        output_format = "delimited"
        output_delimiter="\t"
    elif args.output_format == "csv":
        suffix = ".tsv"
        output_format = "delimited"
        output_delimiter=","
    else:
        raise ValueError(args.output_format)
    if args.output_prefix is None:
        args.output_prefix = os.path.splitext(args.config_file)[0] + ".delineate"
    if args.output_prefix == "-":
        out_name = "<stdout>"
        out = sys.stdout
    else:
        out_name = args.output_prefix + suffix
        out = open(out_name, "w")
    if args.tree_file:
        tree = model.LineageTree.get(
                path=args.tree_file,
                schema=args.tree_format,
                )
    else:
        tree = None
    with out:
        utility.report_configuration(
                config_d,
                tree=tree,
                logger=logger,
                output_file=out,
                output_format = output_format,
                output_delimiter="\t",
                is_fail_on_extra_tree_lineages=not args.ignore_extra_tree_lineages,
                is_fail_on_extra_configuration_lineages=not args.ignore_extra_configuration_lineages,
                )
    logger.info("DELINEATE configuration data written to: '{}'".format(out_name))

if __name__ == "__main__":
    main()

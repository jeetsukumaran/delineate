#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import argparse
import delineate
from delineate import utility
from delineate import model
from delineate import control

"""
Generate normalized JSON/delimited format configuration file from a JSON/delimited data file.
"""

__prog__ = os.path.basename(__file__)
__version__ = "1.0.0"
__description__ = __doc__
__author__ = 'Jeet Sukumaran and Mark T. Holder'
__copyright__ = 'Copyright (C) 2019 Jeet Sukumaran and Mark T. Holder.'

def main():
    parser = argparse.ArgumentParser(description=__doc__,
            parents=[delineate.get_metadata_parser_opts()])
    parser.add_argument("-c", "--constraints-file",
            action="store",
            required=True,
            help="Path to constraints configuration file; can be JSON for columnar tab/comma-delimited file.")
    parser.add_argument("-d", "--config-file-delimiter",
            action="store",
            dest="delimiter",
            default=None,
            help="[Assuming non-JSON columnar file] configuration data delimiting character.")
    parser.add_argument("-t", "--tree_file",
            default=None,
            required=True,
            help="Path to tree file (optional; for validation)")
    parser.add_argument("-f", "--tree-format",
            dest="tree_format",
            type=str,
            default="nexus",
            choices=["nexus", "newick"],
            help="Tree file data format (default='%(default)s').")
    parser.add_argument("--output-format",
            choices=["json", "tsv", "csv"],
            default="tsv",
            help="Format for output file (default='%(default)s').")
    parser.add_argument("-o", "--output-prefix",
            action="store",
            default=None,
            help="Prefix for configuration output files; if provided normalized configuration will be written here.")
    parser.add_argument("--case-sensitive-labels",
            dest="is_case_sensitive",
            action="store_true",
            default=False,
            help="Labels will be treated case-sensitively.")
    parser.add_argument("--ignore-extra-tree-lineages",
            dest="is_fail_on_extra_tree_lineages",
            action="store_false",
            default=True,
            help="Do not complain if not all tree lineages are specified in the configuration file.")
    parser.add_argument("--ignore-extra-configuration-lineages",
            dest="is_fail_on_extra_configuration_lineages",
            action="store_false",
            default=True,
            help="Do not complain if not all configuration lineages are found on the tree.")
    args = parser.parse_args()
    if args.describe:
        delineate.description()
        sys.exit(0)
    if args.version:
        print(delineate.name())
        sys.exit(0)
    if not args.constraints_file:
        sys.exit("Need to specify path to constraints configuration file using '-c' option")
    if args.delimiter is not None and (
            args.delimiter.upper() == "TAB"
            or args.delimiter == '\\t'
            or args.delimiter == '\t'
            ):
        args.delimiter = "\t"
    controller = control.get_controller(
            args=args,
            name="delineate-check",
            logger_kwargs={
                "is_include_name": True,
                "is_include_timestamp": False,
                "log_to_stderr": True,
                "log_to_file": False,
                })
    if not args.output_prefix:
        sys.exit(0)
    # logger = utility.RunLogger(name="delineate-configure",
    #         is_include_name=True,
    #         is_include_timestamp=False,
    #         log_to_stderr=True,
    #         stderr_logging_level=utility.logging.INFO,
    #         log_to_file=False,
    #         file_logging_level=utility.logging.INFO,
    #         )
    # config_d = utility.parse_configuration(
    #         args=args,
    #         delimiter=args.delimiter,
    #         logger=logger)
    if args.output_format == "json":
        suffix = ".json"
        output_format = "json"
        output_delimiter="\t"
    elif args.output_format == "tsv":
        suffix = ".tsv"
        output_format = "delimited"
        output_delimiter="\t"
    elif args.output_format == "csv":
        suffix = ".csv"
        output_format = "delimited"
        output_delimiter=","
    else:
        raise ValueError(args.output_format)
    if args.output_prefix == "-":
        out_desc = " standard output"
        out = sys.stdout
    else:
        out_name = args.output_prefix + suffix
        out = open(out_name, "w")
        out_desc = ": '{}'".format(out_name)
    with out:
        controller.write_configuration(
                output_file=out,
                output_format = output_format,
                output_delimiter="\t",
                )
    if out_desc:
        controller.logger.info("DELINEATE configuration data written to{}".format(out_desc))

if __name__ == "__main__":
    main()

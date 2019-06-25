#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import argparse
from delineate import utility

"""
Generate a JSON format configuration file from a columnar data file.
"""

__prog__ = os.path.basename(__file__)
__version__ = "1.0.0"
__description__ = __doc__
__author__ = 'Jeet Sukumaran and Mark T. Holder'
__copyright__ = 'Copyright (C) 2019 Jeet Sukumaran and Mark T. Holder.'

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("source_filepath",
            action="store",
            help="Source of configuration information.")
    parser.add_argument("-o", "--output-prefix",
            action="store",
            default=None,
            help="Prefix for output file.")
    parser.add_argument("-d", "--delimiter",
            action="store",
            default=None,
            help="Input file delimiter [default=<TABE>].")
    parser.add_argument("--pretty-print",
            action="store_true",
            help="Pretty-print JSON.")
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
    with open(os.path.expandvars(os.path.expanduser(args.source_filepath))) as src:
        config_d = utility.parse_delimited_configuration_file(
                src=src,
                delimiter=args.delimiter,
                logger=logger)
    if args.output_prefix is None:
        out = open(os.path.splitext(args.source_filepath)[0] + ".delineate.json", "w")
    elif args.output_prefix == "-":
        out = sys.stdout
    else:
        out = open(args.output_prefix + ".delineate.json", "w")
    kwargs = {}
    if args.pretty_print:
        kwargs["sort_keys"] = True
        kwargs["indent"] = 4
    with out:
        json.dump(config_d, out, **kwargs)
    if not hasattr(out, "name"):
        out_name = "<stdout>"
    else:
        out_name = out.name
    logger.info("DELINEATE configuration data written to: '{}'".format(out_name))

if __name__ == "__main__":
    main()

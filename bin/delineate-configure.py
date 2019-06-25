#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import csv
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
            default="\t",
            help="Input file delimiter [default=<TABE>].")
    parser.add_argument("--pretty-print",
            action="store_true",
            help="Pretty-print JSON.")
    args = parser.parse_args()
    logger = utility.RunLogger(name="delineate-configure",
            is_include_name=True,
            is_include_timestamp=False)
    with open(os.path.expandvars(os.path.expanduser(args.source_filepath))) as src:
        src_data = csv.DictReader(src,
            delimiter=args.delimiter,
            quoting=csv.QUOTE_NONE,
            )
        fieldname_set = set(src_data.fieldnames)
        found_fields = "\n".join("-    '{}'".format(f) for f in src_data.fieldnames)
        msg = []
        msg.append(("{} fields found in configuration source:".format(len(src_data.fieldnames))))
        for idx, fn in enumerate(src_data.fieldnames):
            msg.append("    [{}/{}] '{}'".format(idx+1, len(src_data.fieldnames), fn))
        logger.info("\n".join(msg))
        for required_field in ("label", "species", "status"):
            if required_field not in fieldname_set:
                logger.error("[delineate-configure] ERROR: Field '{}' not found in configuration source".format(required_field))
                sys.exit(1)
        species_label_map = {}
        known = 0
        unknown = 0
        for entry in src_data:
            if entry["status"] == "1":
                try:
                    species_label_map[entry["species"]].append(entry["label"])
                except KeyError:
                    species_label_map[entry["species"]] = [entry["label"]]
                known += 1
            elif entry["status"] == "0":
                unknown += 1
                pass
            else:
                sys.exit("Unrecognized status: '{}'".format(entry["status"]))
    logger.info("{} lineages in total".format(known + unknown))
    msg = []
    msg.append("{} species defined in total".format(len(species_label_map)))
    for spidx, sp in enumerate(species_label_map):
        msg.append("    [{}/{}] '{}'".format(spidx+1, len(species_label_map), sp))
    msg = "\n".join(msg)
    logger.info(msg)
    logger.info("{} lineages assigned to {} species".format(known, len(species_label_map)))
    logger.info("{} lineages of unknown species affinities".format(unknown))
    species_leafset_constraints = []
    for key in species_label_map:
        species_leafset_constraints.append(species_label_map[key])
    assert len(species_leafset_constraints) == len(species_label_map)
    config_d = {}
    config_d["species_leafset_constraints"] = species_leafset_constraints
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

#! /usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
##
##  Copyright 2018 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
This program does something.
"""

import sys
import os
import re
import argparse
import json
import collections
import scipy.optimize
from delineate import model

__prog__ = os.path.basename(__file__)
__version__ = "1.0.0"
__description__ = __doc__
__author__ = 'Jeet Sukumaran and Mark T. Holder'
__copyright__ = 'Copyright (C) 2018 Jeet Sukumaran and Mark T. Holder.'

def parse_fieldname_and_value(labels):
    if not labels:
        return collections.OrderedDict()
    fieldname_value_map = collections.OrderedDict()
    for label in labels:
        match = re.match(r"\s*(.*?)\s*:\s*(.*)\s*", label)
        if not match:
            raise ValueError("Cannot parse fieldname and label (format required: fieldname:value): {}".format(label))
        fieldname, value = match.groups(0)
        fieldname_value_map[fieldname] = value
    return fieldname_value_map

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__description__)

    source_options = parser.add_argument_group("Source options")

    source_options.add_argument("-t", "--tree-file",
            nargs=1,
            help="Path to tree file.")

    source_options.add_argument("-c", "--config-file",
            help="Path to configuration file.")

    source_options.add_argument("-f", "--format",
            dest="data_format",
            type=str,
            default="nexus",
            choices=["nexus", "newick"],
            help="Input data format (default='%(default)s').")

    output_options = parser.add_argument_group("Output options")
    output_options.add_argument("-i", "--tree-info",
            action="store_true",
            help="Output additional information about the tree",)
    output_options.add_argument("-l", "--label",
            action="append",
            help="Label to append to output (in format <FIELD-NAME>:value;)")
    output_options.add_argument( "--no-header-row",
            action="store_true",
            default=False,
            help="Do not write a header row.")
    output_options.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append to output file if it already exists instead of overwriting.")

    args = parser.parse_args()
    args.separator = "\t"
    tree = model.LineageTree.get(
            path=args.tree_file[0],
            schema=args.data_format,
            )
    with open(args.config_file) as src:
        config = json.load(src)
    initial_speciation_rate = config.pop("initial_speciation_rate", 0.01)
    min_speciation_rate = config.pop("min_speciation_rate", 0.00)
    max_speciation_rate = config.pop("max_speciation_rate", 2.00)
    species_leaf_sets = model._Partition.compile_lookup_key( config["species_leaf_sets"] )
    def f(x, *args):
        tree.speciation_completion_rate = x
        return -1 * tree.calc_joint_probability_of_species(taxon_labels=species_leaf_sets)
    #scipy.optimize.bracket(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0, maxiter=1000)[source]
    brac_res = scipy.optimize.bracket(f,
            xa=min_speciation_rate,
            xb=max_speciation_rate,
            )
    brent_res = scipy.optimize.brent(f, brack=brac_res[:3])
    extra_fields = parse_fieldname_and_value(args.label)
    if not args.no_header_row:
        header_row = []
        if args.tree_info:
            header_row.append("ntips")
            header_row.append("rootAge")
        header_row.extend(extra_fields)
        header_row.append("estSpCompRate")
        sys.stdout.write(args.separator.join(header_row))
        sys.stdout.write("\n")
    row = []
    if args.tree_info:
        row.append("{}".format(len(tree.taxon_namespace)))
        tree.calc_node_ages()
        row.append("{}".format(tree.seed_node.age))
    for field in extra_fields:
        row.append(extra_fields[field])
    row.append("{}".format(brent_res))
    sys.stdout.write(args.separator.join(row))
    sys.stdout.write("\n")

if __name__ == '__main__':
    main()



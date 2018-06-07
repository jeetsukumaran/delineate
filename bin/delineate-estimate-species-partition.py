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
import math
from delineate import model
from delineate import utility

__prog__ = os.path.basename(__file__)
__version__ = "1.0.0"
__description__ = __doc__
__author__ = 'Jeet Sukumaran and Mark T. Holder'
__copyright__ = 'Copyright (C) 2018 Jeet Sukumaran and Mark T. Holder.'

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__description__)
    utility.add_source_options(parser)
    model_options = parser.add_argument_group("Model options")
    model_options.add_argument("-s", "--speciation-completion-rate",
            type=float,
            default=None,
            help="The speciation completion rate (overrides that given in configuration file, if any; default: None).",)
    output_options = parser.add_argument_group("Output options")
    output_options.add_argument("-I", "--tree-info",
            action="store_true",
            help="Output additional information about the tree",)
    utility.add_output_options(parser, output_options=output_options)
    args = parser.parse_args()
    tree = model.LineageTree.get(
            path=args.tree_file,
            schema=args.data_format,
            )
    with open(args.config_file) as src:
        config = json.load(src)
    speciation_completion_rate = config.get("speciation_completion_rate", args.speciation_completion_rate)
    tree.speciation_completion_rate = speciation_completion_rate
    partition_probability_map = tree.calc_label_partition_probability_map()
    for key in partition_probability_map:
        print("{}: {}".format(key, partition_probability_map[key]))
    # species_leaf_sets = model._Partition.compile_lookup_key( config["species_leaf_sets"] )
    # initial_speciation_rate = config.pop("initial_speciation_rate", 0.01)
    # min_speciation_rate = config.pop("min_speciation_rate", 1e-8)
    # max_speciation_rate = config.pop("max_speciation_rate", 2.00)
    # mle = MaximumLikelihoodEstimator(
    #         tree=tree,
    #         species_leaf_sets=species_leaf_sets,
    #         initial_speciation_rate=initial_speciation_rate,
    #         min_speciation_rate=min_speciation_rate,
    #         max_speciation_rate=max_speciation_rate)
    # speciation_completion_rate_estimate, speciation_completion_rate_estimate_lnl = mle.estimate_speciation_rate()
    # extra_fields = utility.parse_fieldname_and_value(args.label)
    # if not args.no_header_row:
    #     header_row = []
    #     if args.tree_info:
    #         header_row.append("numTips")
    #         header_row.append("rootAge")
    #     header_row.extend(extra_fields)
    #     header_row.append("estSpCompRate")
    #     header_row.append("estSpCompRateLnL")
    #     if args.intervals:
    #         header_row.append("ciLow")
    #         header_row.append("ciHigh")
    #     sys.stdout.write(args.separator.join(header_row))
    #     sys.stdout.write("\n")
    # row = []
    # if args.tree_info:
    #     row.append("{}".format(len(tree.taxon_namespace)))
    #     tree.calc_node_ages()
    #     row.append("{}".format(tree.seed_node.age))
    # for field in extra_fields:
    #     row.append(extra_fields[field])
    # row.append("{}".format(speciation_completion_rate_estimate))
    # row.append("{}".format(speciation_completion_rate_estimate_lnl))
    # if args.intervals:
    #     ci_low, ci_high = mle.estimate_confidence_interval(
    #         mle_speciation_rate=speciation_completion_rate_estimate,
    #         max_lnl=speciation_completion_rate_estimate_lnl)
    #     row.append("{}".format(ci_low))
    #     row.append("{}".format(ci_high))
    # sys.stdout.write(args.separator.join(row))
    # sys.stdout.write("\n")

if __name__ == '__main__':
    main()



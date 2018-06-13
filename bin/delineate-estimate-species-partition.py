#! /usr/bin/env python3
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
    args.format = "json"

    tree = model.LineageTree.get(
            path=args.tree_file,
            schema=args.data_format,
            )
    with open(args.config_file) as src:
        config = json.load(src)
    speciation_completion_rate = config.get("speciation_completion_rate", args.speciation_completion_rate)
    if speciation_completion_rate is None:
        raise ValueError("Speciation completion rate must be specified either in configuration file or by command argument '--speciation-completion-rate'")
    tree.speciation_completion_rate = speciation_completion_rate
    if "species_constraints" in config:
        species_constraints = model._Partition.compile_lookup_key( config["species_constraints"] )
        tree.set_node_constraints(species_constraints)
    partition_probability_map = tree.calc_label_partition_probability_map()

    result_dict = collections.OrderedDict()
    # if not args.no_header_row:
    #     header_row = []
    #     if args.tree_info:
    #         header_row.append("numTips")
    #         header_row.append("rootAge")
    #     header_row.extend(extra_fields)
    #     header_row.append("estSpCompRate")
    #     header_row.append("estSpCompRateLnL")
    #     sys.stdout.write(args.separator.join(header_row))
    #     sys.stdout.write("\n")
    row = []
    if args.tree_info:
        result_dict["numTips"] = len(tree.taxon_namespace)
        tree.calc_node_ages()
        result_dict["rootAge"] = "{}".format(tree.seed_node.age)
    extra_fields = utility.parse_fieldname_and_value(args.label)
    if extra_fields:
        result_dict.update(extra_fields)

    species_partition_info = [(k, list(list(s) for s in k), math.log(partition_probability_map[k])) for k in partition_probability_map]
    species_partition_info.sort(key=lambda x: x[2], reverse=True)
    assert species_partition_info[0][2] >= species_partition_info[-1][2]
    result_dict["numPartitions"] = len(species_partition_info)
    max_loglikelihood = species_partition_info[0][2]
    result_dict["lnLMax"] = max_loglikelihood
    result_dict["partitions"] = []
    for key_idx, (key, key_as_list, lnL) in enumerate(species_partition_info):
        p = collections.OrderedDict()
        p["partition"] = key_as_list
        p["lnL"] = lnL
        if lnL + 1.96 >= max_loglikelihood:
            p["isInConfidenceInterval"] = True
        else:
            p["isInConfidenceInterval"] = False
        # sys.stderr.write("{}, {}, {}\n".format(max_loglikelihood, lnL, p["isInConfidenceInterval"]))
        result_dict["partitions"].append(p)
    if args.format == "json-compact":
        json.dump(result_dict, sys.stdout)
    else:
        json.dump(result_dict, sys.stdout, indent=4, separators=(',', ': '))
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



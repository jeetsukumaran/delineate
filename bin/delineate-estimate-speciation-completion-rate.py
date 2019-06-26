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
from delineate import model
from delineate import estimate
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

    estimation_options = parser.add_argument_group("Estimation options")
    estimation_options.add_argument("-i", "--intervals", "--confidence-intervals",
            action="store_true",
            help="Calculate confidence intervals.",)
    estimation_options.add_argument("-u", "--underflow-protection",
            action="store_true",
            default=False,
            help="Try to protect against underflow by using special number handling classes (slow).",)

    output_options = parser.add_argument_group("Output options")
    output_options.add_argument("-I", "--tree-info",
            action="store_true",
            help="Output additional information about the tree",)
    utility.add_output_options(parser, output_options=output_options)

    args = parser.parse_args()
    logger = utility.RunLogger(name="delineate-estimate")
    args.output_field_separator = "\t"
    tree = model.LineageTree.get(
            path=args.tree_file,
            schema=args.tree_format,
            )
    tree.is_use_decimal_value_type = args.underflow_protection
    config = utility.parse_configuration(
            args=args,
            logger=logger)
    utility.report_configuration(
            config_d=config,
            tree=tree,
            logger=logger,
            output_file=None,
            is_fail_on_extra_tree_lineages=False,
            is_fail_on_extra_configuration_lineages=True,
            )
    if "species_leafsets" in config and utility.SPECIES_LEAFSET_CONSTRAINTS_KEY in config:
        sys.exit("Both 'species_leafsets' and '{}' specified in configuration".format(utility.SPECIES_LEAFSET_CONSTRAINTS_KEY))
    elif "species_leafsets" in config:
        key = "species_leafsets"
    elif utility.SPECIES_LEAFSET_CONSTRAINTS_KEY in config:
        key = utility.SPECIES_LEAFSET_CONSTRAINTS_KEY
    else:
        sys.exit("Neither 'species_leafsets' nor '{}' specified in configuration".format(utility.SPECIES_LEAFSET_CONSTRAINTS_KEY))
    species_leafset_labels = model._Partition.compile_lookup_key( config[key] )
    initial_speciation_rate = config.pop("initial_speciation_rate", 0.01)
    min_speciation_rate = config.pop("min_speciation_rate", 1e-8)
    max_speciation_rate = config.pop("max_speciation_rate", 2.00)
    mle = estimate.SpeciationCompletionRateMaximumLikelihoodEstimator(
            tree=tree,
            species_leafset_labels=species_leafset_labels,
            initial_speciation_rate=initial_speciation_rate,
            min_speciation_rate=min_speciation_rate,
            max_speciation_rate=max_speciation_rate)
    speciation_completion_rate_estimate, speciation_completion_rate_estimate_lnl = mle.estimate_speciation_rate()
    extra_fields = utility.parse_fieldname_and_value(args.label)
    if not args.no_header_row:
        header_row = []
        if args.tree_info:
            header_row.append("num_tips")
            header_row.append("root_age")
        header_row.extend(extra_fields)
        header_row.append("speciation_completion_rate")
        header_row.append("speciation_completion_rate_estimate_lnl")
        if args.intervals:
            header_row.append("ci_low")
            header_row.append("ci_high")
        sys.stdout.write(args.output_field_separator.join(header_row))
        sys.stdout.write("\n")
    row = []
    if args.tree_info:
        row.append("{}".format(len(tree.taxon_namespace)))
        tree.calc_node_ages()
        row.append("{}".format(tree.seed_node.age))
    for field in extra_fields:
        row.append(extra_fields[field])
    row.append("{}".format(speciation_completion_rate_estimate))
    row.append("{}".format(speciation_completion_rate_estimate_lnl))
    if args.intervals:
        ci_low, ci_high = mle.estimate_confidence_interval(
            mle_speciation_rate=speciation_completion_rate_estimate,
            max_lnl=speciation_completion_rate_estimate_lnl)
        row.append("{}".format(ci_low))
        row.append("{}".format(ci_high))
    sys.stdout.write(args.output_field_separator.join(row))
    sys.stdout.write("\n")

if __name__ == '__main__':
    main()



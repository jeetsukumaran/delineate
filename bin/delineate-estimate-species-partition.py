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
Calculate probabilities of partitions of population tree leafsets under a speciation model.
"""

import sys
import os
import re
import argparse
import json
import collections
import math
import itertools
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
    estimation_options.add_argument("-u", "--underflow-protection",
            action="store_true",
            default=False,
            help="Try to protect against underflow by using special number handling classes (slow).",)
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
    args.output_format = "json-compact"
    logger = utility.RunLogger(name="delineate-estimate")
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

    speciation_completion_rate = config.get("speciation_completion_rate", args.speciation_completion_rate)
    species_leafset_constraint_labels = config.get(utility.SPECIES_LEAFSET_CONSTRAINTS_KEY, None)
    if species_leafset_constraint_labels is not None:
        species_constraints = model._Partition.compile_lookup_key(species_leafset_constraint_labels)
        tree.set_node_constraints(species_constraints)

    if speciation_completion_rate is None:
        if species_leafset_constraint_labels is not None:
            constraint_tip_labels = list(itertools.chain(*species_leafset_constraint_labels))
            induced_tree = tree.extract_tree_with_taxa_labels(
                    labels=constraint_tip_labels)
            mle = estimate.SpeciationCompletionRateMaximumLikelihoodEstimator(
                    tree=induced_tree,
                    species_leafset_labels=species_constraints,
                    initial_speciation_rate=config.get("initial_speciation_rate", 0.01),
                    min_speciation_rate=config.get("min_speciation_rate", 1e-8),
                    max_speciation_rate=config.get("max_speciation_rate", 2.00))
            speciation_completion_rate_estimate, speciation_completion_rate_estimate_lnl = mle.estimate_speciation_rate()
            speciation_completion_rate = speciation_completion_rate_estimate
            speciation_completion_rate_source = "estimated"
    else:
        speciation_completion_rate_source = "specified"
        speciation_completion_rate_estimate_lnl = 0.0
    if speciation_completion_rate is None:
        raise ValueError("Speciation completion rate must be specified either in configuration file or by command argument '--speciation-completion-rate'")
    tree.speciation_completion_rate = speciation_completion_rate

    partition_probability_map = tree.calc_label_partition_probability_map()

    result_dict = collections.OrderedDict()
    row = []
    if args.tree_info:
        result_dict["num_tips"] = len(tree.taxon_namespace)
        tree.calc_node_ages()
        result_dict["root_age"] = "{}".format(tree.seed_node.age)
    extra_fields = utility.parse_fieldname_and_value(args.label)
    if extra_fields:
        result_dict.update(extra_fields)
    result_dict["speciation_completion_rate"] = speciation_completion_rate
    result_dict["speciation_completion_rate_source"] = speciation_completion_rate_source
    result_dict["speciation_completion_rate_estimate_lnl"] = speciation_completion_rate_estimate_lnl

    species_partition_info = []
    probability_index = 2
    for k in partition_probability_map:
        kentry = (
                k,
                list(list(s) for s in k),
                partition_probability_map[k],
                )
        species_partition_info.append(kentry)
    # species_partition_info = [(
    #         k,
    #         list(list(s) for s in k),
    #         partition_probability_map[k],
    #         math.log(partition_probability_map[k])) for k in partition_probability_map]
    species_partition_info.sort(key=lambda x: x[probability_index], reverse=True)
    assert species_partition_info[0][probability_index] >= species_partition_info[-1][probability_index], \
            sys.stderr.write("{} >= {}: False\n".format(species_partition_info[0][probability_index], species_partition_info[-1][probability_index]))
    result_dict["num_partitions"] = len(species_partition_info)
    # max_log_likelihood = species_partition_info[0][log_likelihood_index]
    # result_dict["max_log_likelihood"] = max_log_likelihood
    result_dict["num_partitions_in_confidence_interval"] = 0
    result_dict["partitions"] = []
    cumulative_probability = tree.as_working_value_type(0.0)
    cumulative_probability_given_constr = tree.as_working_value_type(0.0)
    num_partitions_in_confidence_interval = 0
    cond_prob = sum([i[probability_index] for i in species_partition_info])
    if cond_prob == 0:
        raise ValueError("0 conditional probability")
    try:
        ln_cond_prob = math.log(cond_prob)
    except ValueError:
        ln_cond_prob = float("nan")
    for key_idx, (key, key_as_list, prob) in enumerate(species_partition_info):
        p = collections.OrderedDict()
        p["species_leafsets"] = key_as_list
        # try:
        #     p["log_probability"] = math.log(prob)
        # except ValueError:
        #     p["log_probability"] = float("nan")
        p["probability"] = tree.as_float(prob)
        bpc = prob / cond_prob
        p["probability_given_constraints"] = tree.as_float(bpc)
        # p["log_probability_given_constraints"] = tree.as_float(lnL - ln_cond_prob)

        # need to check this before summing cumulative probability, otherwise
        # any single partition with probability > 0.95 will incorrectly be
        # flagged as being outside the confidence interval.
        if cumulative_probability_given_constr <= 0.95 and bpc > 0.0:
            p["is_in_confidence_interval"] = True
            num_partitions_in_confidence_interval += 1
        else:
            p["is_in_confidence_interval"] = False

        cumulative_probability += prob
        cumulative_probability_given_constr += bpc
        p["cumulative_probability"] = tree.as_float(cumulative_probability)
        p["cumulative_probability_given_constraints"] = tree.as_float(cumulative_probability_given_constr)

        result_dict["partitions"].append(p)
    result_dict["num_partitions_in_confidence_interval"] = num_partitions_in_confidence_interval
    if args.output_format == "json-compact":
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



#! /usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
##
##  Copyright 2019 Jeet Sukumaran and Mark T. Holder.
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
Calculate probabilities of partitions of population tree leafsets into species
under a Protracted Birth Death (PBD) model, optimizing the speciation
completion rate parameter as required.
"""

# Imports {{{1
import contextlib
import sys
import os
import argparse
import subprocess
import math
import itertools
import collections
import json
from dendropy.dataio import nexusprocessing
from dendropy.model import birthdeath
import delineate
from delineate import model
from delineate import estimate
from delineate import control
from delineate import utility
# }}}1

# Support Functions {{{1

def _node_label_compose_fn(node, key):
    return node.annotations.get_value(key, None)

def _check_labels_in_tree(labels, tree):
    """
    Used to make sure the internal labels are consistent.
    This should be after all normalization of labels is done.
    By this point, all space/underscore issues, casing, etc. should be homogenous or accounted for.
    """
    tree_lineage_labels = set([tx.label for tx in tree.taxon_namespace])
    not_found = []
    for label in labels:
        if label not in tree_lineage_labels:
            not_found.append(label)
    return not_found

def sorted_dict(d):
    d2 = {}
    for k in sorted(d.keys()):
        d2[k] = d[k]
    return d2
# }}}1

# Species Partition Estimation {{{1
def execute_species_partition_estimation(args):
    controller = control.get_controller(
            name="delineate-estimate",
            args=args)
    if args.speciation_completion_rate is not None:
        controller.speciation_completion_rate = args.speciation_completion_rate
    speciation_completion_rate = controller.speciation_completion_rate
    tree = controller.tree
    tree.birth_rate = birthdeath.fit_pure_birth_model_to_tree(tree=tree)["birth_rate"]
    if controller.has_species_constraints:
        species_constraints = model._Partition.compile_lookup_key(controller.species_leafset_constraint_labels)
        tree.set_node_constraints(species_constraints)
        check_labels = _check_labels_in_tree(
                labels=controller.constrained_lineage_leaf_labels,
                tree=tree)
        if check_labels:
            raise ValueError("Lineage labels not normalized or invalid: {}".format(check_labels))
        induced_tree = tree.extract_tree_with_taxa_labels(
                labels=controller.constrained_lineage_leaf_labels)
        speciation_completion_rate_estimation_initial = controller.speciation_completion_rate_estimation_initial
        speciation_completion_rate_estimation_min = controller.speciation_completion_rate_estimation_min
        speciation_completion_rate_estimation_max = controller.speciation_completion_rate_estimation_max
        mle = estimate.SpeciationCompletionRateMaximumLikelihoodEstimator(
                tree=induced_tree,
                species_leafset_labels=species_constraints,
                initial_speciation_rate=speciation_completion_rate_estimation_initial,
                min_speciation_rate=speciation_completion_rate_estimation_min,
                max_speciation_rate=speciation_completion_rate_estimation_max,
                )
        constrained_tree = induced_tree.as_string(
                schema="newick",
                )
        speciation_completion_rate_estimate, speciation_completion_rate_estimate_lnl = mle.estimate_speciation_rate()
        speciation_completion_rate = speciation_completion_rate_estimate
        speciation_completion_rate_source = "estimated"
        # species_constraints_desc = controller.species_leafset_constraint_labels
        try:
            species_constraints_desc = dict(controller.config_d["configuration_table"])
            try:
                del species_constraints_desc["lineages"]
            except KeyError:
                pass
        except KeyError:
            # probably given a JSON file here;
            # TODO: pull out what we can
            # species_constrained_lineage_map
            # constrained_lineage_species_map
            species_constrained_lineage_map = {}
            constrained_lineage_species_map = {}
            slcd = controller.config_d["species_leafset_constraints"]
            for spidx in range(len(slcd)):
                species_name = "ConstrainedSp{:03d}".format(spidx+1)
                species_constrained_lineage_map[species_name] = list(slcd[spidx])
                for lineage in slcd[spidx]:
                    constrained_lineage_species_map[lineage] = species_name

            species_constraints_desc = {
                    "species_constrained_lineage_map": species_constrained_lineage_map,
                    "constrained_lineage_species_map": constrained_lineage_species_map,
            }
    else:
        speciation_completion_rate_estimation_initial = None
        speciation_completion_rate_estimation_min = None
        speciation_completion_rate_estimation_max = None
        constrained_tree = None
        species_constraints_desc = None
    if speciation_completion_rate is None and not controller.has_species_constraints:
        controller.logger.error("ERROR: With no constraints given, speciation completion rate must be specified either by command argument '--speciation-completion-rate'")
        controller.logger.critical("Terminating due to error")
        sys.exit(1)
    elif speciation_completion_rate is not None and not controller.has_species_constraints:
        speciation_completion_rate_source = "specified"
        speciation_completion_rate_estimate_lnl = 0.0
    tree.speciation_completion_rate = speciation_completion_rate

    partition_probability_map = tree.calc_label_partition_probability_map()

    result_dict = collections.OrderedDict()
    row = []
    if args.tree_info:
        result_dict["num_tips"] = len(tree.taxon_namespace)
        tree.calc_node_ages()
        result_dict["root_age"] = "{}".format(tree.seed_node.age)
    extra_fields = utility.parse_fieldname_and_value(args.extra_info_field_value)
    if extra_fields:
        result_dict.update(extra_fields)
    result_dict["lineages"] = [txn.label for txn in tree.taxon_namespace]
    result_dict["speciation_completion_rate"] = speciation_completion_rate
    result_dict["speciation_completion_rate_source"] = speciation_completion_rate_source
    result_dict["speciation_completion_rate_estimate_lnl"] = speciation_completion_rate_estimate_lnl
    result_dict["speciation_completion_rate_estimation_initial"] = speciation_completion_rate_estimation_initial
    result_dict["speciation_completion_rate_estimation_min"] = speciation_completion_rate_estimation_min
    result_dict["speciation_completion_rate_estimation_max"] = speciation_completion_rate_estimation_max
    result_dict["lineage_tree_birth_rate"] = controller.tree.birth_rate
    result_dict["constrained_lineage_tree"] = constrained_tree
    result_dict["species_constraints"] = species_constraints_desc
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
    result_dict["report_mle_only"] = args.report_mle_only
    result_dict["report_constrained_probability_threshold"] = args.report_constrained_probability_threshold
    result_dict["report_constrained_cumulative_probability_threshold"] = args.report_constrained_cumulative_probability_threshold
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

    dest_tree_files = []
    if not args.no_primary_tree_file:
        treefile_path, treefile = open_output_file(
                args=args,
                suffix=".delimitation-results",
                extension="trees",)
        dest_tree_files.append(treefile)
    if args.truncated_tree_file_size > 0:
        trunc_treefile_path, trunc_treefile = open_output_file(
                args=args,
                suffix=".delimitation-results",
                extension="trunc.trees",)
        trunc_treefile_idx = len(dest_tree_files)
        dest_tree_files.append(trunc_treefile)
    else:
        trunc_treefile_idx = None
    leaf_nodes = [nd for nd in tree.leaf_node_iter()]
    if args.store_relabeled_trees is not None:
        relabeled_trees_path, rltf = open_output_file(
            args=args,
            suffix=".delimitation-results.relabeled",
            extension="trees",)
        if args.store_relabeled_trees == "species-lineage":
            relabel_key = "species-lineage"
        elif args.store_relabeled_trees == "lineage-species":
            relabel_key = "lineage-species"
        else:
            raise ValueError(args.store_relabeled_trees)

    for d in dest_tree_files:
        d.is_active = True
    def _wt(s):
        for d in dest_tree_files:
            if d.is_active:
                d.write(s)
    # with ExitStack() as stack:
    #     files = [stack.enter_context(open(fname)) for fname in filenames]
    #     Do something with "files"
    with contextlib.ExitStack() as stack:
        for dest_tree_file in dest_tree_files:
            stack.enter_context(dest_tree_file)
            tree.taxon_namespace.write(file=dest_tree_file, schema="nexus")
        _wt("BEGIN TREES;\n")
        taxon_token_map = {}
        translate_statement = []
        for taxon_idx, taxon in enumerate(tree.taxon_namespace):
            taxon_token_map[taxon] = str(taxon_idx + 1)
            label = nexusprocessing.escape_nexus_token(str(taxon.label),
                    preserve_spaces=True,
                    quote_underscores=True,
                    )
            translate_statement.append("             {} {}".format(taxon_token_map[taxon], label))
        if not args.no_translate_tree_tokens:
            _wt("        Translate\n")
            _wt("{}\n             ;\n".format(",\n".join(translate_statement)))
        for key_idx, (key, key_as_list, prob) in enumerate(species_partition_info):
            if args.report_mle_only and key_idx > 0:
                break
            if key_idx == args.truncated_tree_file_size:
                dest_tree_files[trunc_treefile_idx].is_active = False
            p = collections.OrderedDict()
            lineage_species_name_map = controller.compile_postanalysis_lineage_species_name_map(postanalysis_species_leafset_labels=key_as_list)
            p["lineage_species_name_map"] = sorted_dict(lineage_species_name_map)
            species_lineage_name_map = {}
            for lineage in lineage_species_name_map:
                try:
                    species_lineage_name_map[lineage_species_name_map[lineage]].append(lineage)
                except KeyError:
                    species_lineage_name_map[lineage_species_name_map[lineage]] = [lineage]
            p["species_leafsets"] = sorted_dict(species_lineage_name_map)
            p["constrained_probability"] = None
            p["constrained_cumulative_probability"] = None
            p["is_in_confidence_interval"] = None
            p["unconstrained_probability"] = None
            p["unconstrained_cumulative_probability"] = None
            # try:
            #     p["log_probability"] = math.log(prob)
            # except ValueError:
            #     p["log_probability"] = float("nan")
            p["unconstrained_probability"] = tree.as_float(prob)
            bpc = prob / cond_prob
            p["constrained_probability"] = tree.as_float(bpc)
            if args.report_constrained_probability_threshold is not None and p["constrained_probability"] < args.report_constrained_probability_threshold:
                break
            # if args.report_unconstrained_probability_threshold is not None and p["unconstrained_probability"] < args.report_unconstrained_probability_threshold:
            #     break
            # p["log_probability_given_constraints"] = tree.as_float(lnL - ln_cond_prob)

            # need to check this before summing cumulative probability, otherwise
            # any single partition with probability > 0.95 will incorrectly be
            # flagged as being outside the confidence interval.
            if cumulative_probability_given_constr <= 0.95 and bpc > 0.0:
                p["is_in_confidence_interval"] = True
                num_partitions_in_confidence_interval += 1
            else:
                p["is_in_confidence_interval"] = False
            if args.report_constrained_cumulative_probability_threshold is not None and cumulative_probability_given_constr > args.report_constrained_cumulative_probability_threshold:
                break

            cumulative_probability += prob
            cumulative_probability_given_constr += bpc
            p["unconstrained_cumulative_probability"] = tree.as_float(cumulative_probability)
            p["constrained_cumulative_probability"] = tree.as_float(cumulative_probability_given_constr)

            for key in (
                    "constrained_probability",
                    "constrained_cumulative_probability",
                    "unconstrained_probability",
                    "unconstrained_cumulative_probability",
                    "is_in_confidence_interval",
                    ):
                tree.annotations[key] = p[key]
            for nd in leaf_nodes:
                nd.annotations["species"] = lineage_species_name_map[nd.taxon.label]
                if args.figtree_display_label != "clear":
                    nd.annotations["lineage-species"] = "{} ({})".format(nd.taxon.label, lineage_species_name_map[nd.taxon.label])
                    nd.annotations["species-lineage"] = "{} ({})".format(lineage_species_name_map[nd.taxon.label], nd.taxon.label)
                if nd.taxon.label in controller.registry.preanalysis_constrained_lineage_species_map:
                    nd.annotations["status"] = "constrained"
                else:
                    nd.annotations["status"] = "inferred"
            _wt("    Tree {} = {}".format(
                key_idx + 1,
                tree.as_string(
                    schema="newick",
                    suppress_annotations=False,
                    taxon_token_map=taxon_token_map,
                    )))
            if args.store_relabeled_trees is not None:
                rltf.write(tree.as_string(
                    schema="newick",
                    suppress_annotations=True,
                    node_label_compose_fn=lambda node: _node_label_compose_fn(node, relabel_key),
                    ))
            result_dict["partitions"].append(p)
            # print("--")
            # print(lineage_species_name_map)
        dest_tree_files[trunc_treefile_idx].is_active = True
        _wt("END TREES;\n")
        if args.figtree_display_label != "clear":
            _wt("\n\nbegin figtree;\n")
            if args.figtree_display_label == "species":
                _wt('set tipLabels.displayAttribute="species";\n');
            elif args.figtree_display_label == "species-lineage":
                _wt('set tipLabels.displayAttribute="species-lineage";\n');
            elif args.figtree_display_label == "lineage-species":
                _wt('set tipLabels.displayAttribute="lineage-species";\n');
            _wt("end figtree;\n")

    result_dict["num_partitions_in_confidence_interval"] = num_partitions_in_confidence_interval
    if True: # args.write_json
        out_path, outf = open_output_file(
                args=args,
                suffix=".delimitation-results",
                extension="json")
        with outf:
            if False: #args.output_format == "json-compact":
                json.dump(result_dict, outf)
            else:
                json.dump(result_dict, outf, indent=4, separators=(',', ': '))
            controller.logger.info("JSON-formatted results written to: '{}'".format(out_path))
    controller.logger.info("Operation complete")
    controller.logger.info("Terminating normally")
# }}}1

# Speciation Completion Rate Parameter Estimation {{{1
def execute_speciation_completion_rate_estimation(args):
    controller = control.get_controller(
            name="delineate-estimate",
            args=args)
    tree = controller.tree
    if not controller.has_species_constraints:
        controller.logger.error("ERROR: Species constraints must be fully specifed for rates estimation, but no constraints were found.")
        controller.logger.critical("Terminating due to error")
        sys.exit(1)
    # species_constraints = model._Partition.compile_lookup_key(controller.species_leafset_constraint_labels)
    species_leafset_labels = model._Partition.compile_lookup_key(controller.species_leafset_constraint_labels)
    speciation_completion_rate_estimation_initial = controller.speciation_completion_rate_estimation_initial
    speciation_completion_rate_estimation_min = controller.speciation_completion_rate_estimation_min
    speciation_completion_rate_estimation_max = controller.speciation_completion_rate_estimation_max
    controller.logger.info("Speciation completion rate estimation window minimum: {}".format(speciation_completion_rate_estimation_min))
    controller.logger.info("Speciation completion rate estimation window maximum: {}".format(speciation_completion_rate_estimation_max))
    controller.logger.info("Speciation completion rate estimation initial: {}".format(speciation_completion_rate_estimation_initial))
    mle = estimate.SpeciationCompletionRateMaximumLikelihoodEstimator(
            tree=tree,
            species_leafset_labels=species_leafset_labels,
            initial_speciation_rate=speciation_completion_rate_estimation_initial,
            min_speciation_rate=speciation_completion_rate_estimation_min,
            max_speciation_rate=speciation_completion_rate_estimation_max,
            )
    speciation_completion_rate_estimate, speciation_completion_rate_estimate_lnl = mle.estimate_speciation_rate()
    extra_fields = utility.parse_fieldname_and_value(args.extra_info_field_value)
    out_path, out = open_output_file(
            args=args,
            suffix=".rate-results",
            extension="tsv")
    with out:
        output_field_separator = args.output_field_separator
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
            out.write(output_field_separator.join(header_row))
            out.write("\n")
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
        out.write(output_field_separator.join(row))
        out.write("\n")
        controller.logger.info("Results written to: '{}'".format(out_path))
    controller.logger.info("Operation complete")
    controller.logger.info("Terminating normally")
# }}}1

# CLI Support{{{1
def add_source_options(parser, source_options=None):
    if not source_options:
        source_options = parser.add_argument_group("Source Options")
    source_options.add_argument("-t", "--tree-file",
            required=True,
            help="Path to tree file.")
    source_options.add_argument("-c", "--constraints",
            dest="constraints_file",
            help="Path to constraints file.")
    source_options.add_argument("-f", "--tree-format",
            dest="tree_format",
            type=str,
            default="nexus",
            choices=["nexus", "newick"],
            help="Tree file data format (default='%(default)s').")
    parser.add_argument(
            "--underscores-to-spaces",
            "--no-preserve-underscores",
            action="store_false",
            dest="preserve_underscores",
            default=True,
            help="Convert underscores to spaces in tree lineage labels (if not protected by quotes). "
                 "in the constraints file will not be modified either way. You should ensure "
                 "consistency in labels between the tree file and the constraints file."
            )
    parser._source_options = source_options
    return parser

def add_output_options(parser, output_options=None):
    if not output_options:
        output_options = parser.add_argument_group("Output Options")
    output_options.add_argument("-o", "--output-prefix",
            default=None,
            help="Prefix for output file(s).")
    output_options.add_argument("--extra-info",
            dest="extra_info_field_value",
            action="append",
            help="Extra information to append to output (in format "
                 "'<FIELD-NAME>:value'; e.g. '--extra-info taxonomy:accepted')")
    parser._output_options = output_options
    return parser

def add_estimation_options(parser, estimation_options=None):
    if not estimation_options:
        estimation_options = parser.add_argument_group("Estimation Options")
    estimation_options.add_argument("-u", "--underflow-protection",
            action="store_true",
            default=False,
            help="Try to protect against underflow by using special number handling classes (slow).",)
    estimation_options.add_argument(
            "--speciation-completion-rate-estimation-min",
            "--smin",
            metavar="#.##",
            type=float,
            dest="speciation_completion_rate_estimation_min",
            default=1e-8,
            help="If estimating speciation completion rate, minimum boundary for optimization window  [default: %(default)s].")
    estimation_options.add_argument(
            "--speciation-completion-rate-estimation-max",
            "--smax",
            metavar="#.##",
            type=float,
            dest="speciation_completion_rate_estimation_max",
            default=None,
            help="If estimating speciation completion rate, maximum boundary for optimization window [default: 10 x lineage tree pure birth rate].")
    estimation_options.add_argument(
            "--speciation-completion-rate-estimation-initial",
            "--sinit",
            metavar="#.##",
            type=float,
            dest="speciation_completion_rate_estimation_initial",
            default=None,
            help="If estimating speciation completion rate, initial value for optimizer [default: 0.01 x lineage tree pure birth rate].")
    parser._estimation_options = estimation_options
    return estimation_options

def open_output_file(args,
        suffix,
        extension,
        ):
    if args.output_prefix == "-":
        return "<STDOUT>", sys.stdout
    else:
        if args.output_prefix is None:
            args.output_prefix = utility.compose_output_prefix(
                    input_filepath=args.constraints_file,
                    default="delineate",
                    )
        if suffix is None:
            suffix = ""
        outpath = "{}{}.{}".format(
                args.output_prefix,
                suffix,
                extension)
        if getattr(args, "append", False):
            fmode = "a"
        else:
            fmode = "w"
        return outpath, open(outpath, fmode)
# }}}1

# Main CLI {{{1

def main():
    main_parser = argparse.ArgumentParser(description=__doc__,
            parents=[delineate.get_metadata_parser_opts()])
    try:
        main_parser._optionals.title = "Program Options"
    except AttributeError:
        pass
    cmd_parser = main_parser.add_subparsers(title="Commands", dest="command")

    command_cli = collections.OrderedDict()
    command_cli["partitions"] = {
        "parents": [],
        "desc":
                "Given a known population tree and optionally a speciation"
                " completion rate, calculate the probability of different"
                " partitions of population lineages into species, with the"
                " partition of the highest probability corresponding to the"
                " maximum likelihood species delimitation estimate.",
        "func": execute_species_partition_estimation}
    command_cli["rates"] = {
        "parents": [],
        "desc":
                "Given a known population tree and known species partition,"
                " estimate the speciation completion rate",
        "func": execute_speciation_completion_rate_estimation}
    command_parsers = collections.OrderedDict()
    for key in command_cli:
        kwargs = {}
        if command_cli[key]["parents"]:
            kwargs["parents"] = command_cli[key]["parents"]
        kwargs["help"] = command_cli[key]["desc"]
        kwargs["description"] = kwargs["help"]
        cmd_options_parser = cmd_parser.add_parser(
            key,
            **kwargs)
        cmd_options_parser.set_defaults(func=command_cli[key]["func"])
        try:
            cmd_options_parser._optionals.title = "Command Options"
            cmd_options_parser._positionals.title = "Command Targets"
        except AttributeError:
            pass
        add_source_options(cmd_options_parser)
        add_estimation_options(cmd_options_parser)
        add_output_options(cmd_options_parser)
        command_parsers[key] = cmd_options_parser

    ## setup partition options
    c1_parser = command_parsers["partitions"]
    c1_parser._estimation_options.add_argument("-s", "--speciation-completion-rate",
            type=float,
            default=None,
            help="The speciation completion rate. If specified, then the speciation completion rate "
                 " will *not* be estimated, but fixed to this.",)

    report_options = c1_parser.add_argument_group("Report Options")
    report_options.add_argument("-P", "--report-cumulative-probability-threshold",
            dest="report_constrained_cumulative_probability_threshold",
            metavar="#.##",
            default=None,
            type=float,
            help="Do not report on partitions outside of this constrained (conditional) cumulative probability.")
    report_options.add_argument("--report-mle-only",
            dest="report_mle_only",
            action="store_true",
            default=False,
            help="Only report maximum likelihood estimate.",)
    report_options.add_argument("-p", "--report-probability-threshold",
            dest="report_constrained_probability_threshold",
            metavar="#.##",
            default=None,
            type=float,
            help="Do not report on partitions with individual constrained (conditional) probability below this threshold.")

    output_options = c1_parser._output_options
    output_options.add_argument("-I", "--tree-info",
            action="store_true",
            help="Output additional information about the tree.",)
    trees_options = c1_parser.add_argument_group("Tree Output Options")
    trees_options.add_argument("--no-translate-tree-tokens",
            action="store_true",
            default=False,
            help="Write tree statements using full taxon names rather than numerical indices.")
    trees_options.add_argument("-l", "--figtree-display-label",
            choices=["lineage", "species", "lineage-species", "species-lineage", "clear"],
            default="species-lineage",
            help="Default label to display in FigTree for the trees in the main result file.",)
    trees_options.add_argument("--store-relabeled-trees",
            default=None,
            choices=["lineage-species", "species-lineage"],
            help="Create an additional result file, where trees are written with their labels actually changed to the option selected here.")
    trees_options.add_argument("--truncated-tree-file-size",
            default=10,
            type=int,
            help="Maximum number of trees to include in truncated tree file (for quick reference with larger datasets); default=%(default)s).")
    trees_options.add_argument("--no-primary-tree-file",
            action="store_true",
            default=False,
            help="Do not save primary tree file (truncated tree file will still be written if size > 0).")

    # output_options.add_argument("--output-format",
    #         default="json-compact",
    #         choices=["json", "json-compact"],
    #         help="Format for output file.",)

    ## setup rate estimation options
    c2_parser = command_parsers["rates"]
    estimation_options = c2_parser._estimation_options
    estimation_options.add_argument("-i", "--intervals", "--confidence-intervals",
            action="store_true",
            help="Calculate confidence intervals.",)
    output_options = c2_parser._output_options
    output_options.add_argument( "--no-header-row",
            action="store_true",
            default=False,
            help="Do not write a header row.")
    output_options.add_argument( "--field-separator",
            default="\t",
            dest="output_field_separator",
            help="Field separator or delimiter character [default: tab].")
    output_options.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append to output file if it already exists instead of overwriting.")
    output_options.add_argument("-I", "--tree-info",
            action="store_true",
            help="Output additional information about the tree",)

    if len(sys.argv) == 1:
        msg = []
        msg.append("Please specify one of the following operations:")
        msg.append("")
        msg.append("    delineate-estimate partitions")
        msg.append("    delineate-estimate rates")
        msg.append("")
        msg.append("Or type '--help' for help:")
        msg.append("")
        msg.append("    delineate-estimate --help")
        msg.append("    delineate-estimate partitions --help")
        msg.append("    delineate-estimate rates --help")
        msg.append("")
        msg = "\n".join(msg)
        sys.exit(msg)

    args = main_parser.parse_args()
    if args.describe:
        delineate.description()
        sys.exit(0)
    if args.version:
        print(delineate.name())
        sys.exit(0)
    # controller.is_quiet = getattr(args, "quiet", True)

    if hasattr(args, 'func'):
        args.func(args=args)
    else:
        pass
        # args.status_targets = "short"
        # execute_status_command(controller=controller, args=args)
    sys.exit(0)

if __name__ == '__main__':
    main()

# }}}1

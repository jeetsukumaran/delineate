#! /usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
##
##  Copyright 2020 Jeet Sukumaran and Mark T. Holder.
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
Summarize the results of a BP&P "A10" analysis (delimitation with a
guide tree). Produces the labels of populations for each individual
sequence under different posterior probabilities, as well as a BEAST
"traits" file for each of those probabilities allowing for a StarBeast2
analysis.
"""

usage_examples = [
    "Usage examples:",
    "",
    " $ delineate-bppsum -i bpprun*.imap -r bpprun*.out.txt"
    "",
    " $ delineate-bppsum --imap \\",
    "                           subtree1/st1.imap.txt \\",
    "                           subtree2/st2.imap.txt \\",
    "                           subtree2/st3.imap.txt \\",
    "                   --results \\",
    "                           subtree1/st1.out.txt \\",
    "                           subtree2/st2.out.txt \\",
    "                           subtree3/st3.out.txt \\",
]

import sys
import os
import argparse
import pathlib
import itertools
import dendropy
import csv
from delineate import utility

def flatten(s):
    return itertools.chain.from_iterable(s)

def compile_filepaths(files):
    files = sorted(set(os.path.expanduser(os.path.expandvars(f)) for f in flatten(files)))
    return files

class LineageEntry(object):

    def __init__(self,
            label,
            candidate_population_label,
            imap_filepath):
        self.label = label
        self.candidate_population_label = candidate_population_label
        self.imap_filepath = imap_filepath

class BppSummarizer(object):

    def __init__(self,
            population_label_prefix,
            stderr_logging_level):
        self.logger = utility.RunLogger(
            name="delineate-bppsum",
            is_include_name=True,
            is_include_timestamp=False,
            is_log_to_stderr=True,
            stderr_logging_level=stderr_logging_level,
            is_log_to_file=False,
            file_logging_level=utility.logging.INFO,
            )
        self.population_label_prefix = population_label_prefix
        self.lineages = {}
        self.candidate_population_map = {}
        self.trees = dendropy.TreeList()

    def read_imap(self, files):
        src_files = compile_filepaths(files)
        self.logger.info("{} BPP 'imap' files specified".format(len(src_files)))
        for src_idx, src_filepath in enumerate(src_files):
            self.logger.info("- Reading mapping file {:3d} of {}: {}".format(src_idx+1, len(src_files), src_filepath))
            self._process_imap_data(imap_filepath=src_filepath)

    def _process_imap_data(self, imap_filepath):
        try:
            src = open(imap_filepath)
        except IOError:
            self.logger.error("ERROR opening BPP IMAP file: '{}'".format(imap_filepath))
            sys.exit(1)
        num_lineages = 0
        num_pops = 0
        with src:
            for line in src:
                if not line:
                    continue
                parts = line.split()
                lineage_label = parts[0]
                population = parts[1]
                num_lineages += 1
                if lineage_label in self.lineages:
                    sys.exit("Lineage '{}' already registered from file: '{}'\n".format(
                        lineage_label,
                        self.lineages[lineage_label].imap_filepath))
                self.lineages[lineage_label] = LineageEntry(
                        label=lineage_label,
                        candidate_population_label=population,
                        imap_filepath=imap_filepath,
                        )
                try:
                    self.candidate_population_map[population].append(lineage_label)
                except KeyError:
                    self.candidate_population_map[population] = [lineage_label]
                    num_pops += 1
        self.logger.info("  - ({} lineages, {} candidate populations)".format(
            num_lineages,
            num_pops))

    def read_bpp_output(self, files):
        src_files = compile_filepaths(files)
        self.logger.info("{} BPP output files specified".format(len(src_files)))
        for src_idx, src_filepath in enumerate(src_files):
            self.logger.info("- Reading output file {:3d} of {}: {}".format(src_idx+1, len(src_files), src_filepath))
            self._process_bpp_output_data(out_filepath=src_filepath)
        tree_lineages = set(t.label for t in self.trees.taxon_namespace)
        imap_pops = set(self.candidate_population_map.keys())
        if tree_lineages != imap_pops:
            s = imap_pops - tree_lineages
            msg = ["ERROR: Following {} candidate populations not found in any BPP results tree".format(len(s))]
            for idx, lineage in enumerate(s):
                msg.append("                   - {:3d}: '{}' (defined in '{}')".format(
                    idx+1,
                    lineage,
                    self.lineages[self.candidate_population_map[lineage][0]].imap_filepath))
            self.logger.error("\n".join(msg))
            sys.exit(1)

    def _process_bpp_output_data(self, out_filepath):
        try:
            src = open(out_filepath)
        except IOError:
            self.logger.error("ERROR opening BPP output file: '{}'".format(out_filepath))
            sys.exit(1)
        with src:
            tree_str = src.readlines()[-1].replace("#","")
        tree0 = dendropy.Tree.get(
                data=tree_str,
                schema="newick",
                rooting="force-rooted",
                suppress_external_node_taxa=False,
                suppress_internal_node_taxa=True,
                preserve_underscores=True,
                taxon_namespace=self.trees.taxon_namespace,
                )
        self.trees.append(tree0)
        tree0.encode_bipartitions()
        tree0.src_filepath = out_filepath
        num_lineages = len([nd for nd in tree0.leaf_node_iter()])
        for nd in tree0.preorder_node_iter():
            if getattr(nd, "is_skip_subtree", False):
                continue
            if nd.is_leaf():
                nd.pp = 1.0
                if nd.taxon.label not in self.candidate_population_map:
                    self.logger.error("ERROR: Candidate population lineage '{}' (found in '{}') not defined in any BPP imap file".format(
                        nd.taxon.label,
                        out_filepath))
                    sys.exit(1)
            elif nd.label:
                if nd.label.startswith("#"):
                    nd.pp = float(nd.label[1:])
                else:
                    nd.pp = float(nd.label)
            else:
                nd.pp = 0.0
        self.logger.info("  - ({} candidate populations)".format(num_lineages))

    def calc_populations(self,
            output_prefix,
            population_probability_thresholds,
            field_delimiter,):
        pp_lineage_maps = {}
        for population_probability_threshold in sorted(set(population_probability_thresholds)):
            lineage_population_map = self.calc_lineage_population_map(population_probability_threshold=population_probability_threshold)
            pp_key = "p{:03d}".format(int(population_probability_threshold*100))
            pp_lineage_maps[pp_key] = lineage_population_map
            with open(output_prefix + ".sb2-traits.p{:03d}.txt".format(int(population_probability_threshold*100)), "w") as out:
                out.write("trait\tspecies\n")
                for lineage in sorted(lineage_population_map.keys()):
                    out.write("{:<20}{}{}\n".format(lineage, field_delimiter, lineage_population_map[lineage]))
        with open(output_prefix + ".summary.csv", "w") as out:
            ppkeys = sorted(pp_lineage_maps.keys())
            header_row = ["sample"] + ppkeys
            out.write(field_delimiter.join(header_row))
            out.write("\n")
            for lineage in sorted(self.lineages.keys()):
                row = [lineage]
                for pp in ppkeys:
                    row.append(pp_lineage_maps[pp][lineage])
                out.write(field_delimiter.join(row))
                out.write("\n")

    def calc_lineage_population_map(self, population_probability_threshold):
        new_pop_idx = 1
        lineage_population_map = {}
        for tree in self.trees:
            self.logger.debug("Collapsing populations with <{} probability on tree from: '{}'".format(
                population_probability_threshold,
                tree.src_filepath))
            tree0 = dendropy.Tree(tree)
            tree0.encode_bipartitions()
            for nd in tree0.preorder_node_iter():
                if getattr(nd, "is_skip_subtree", False):
                    continue
                nd.pp = tree.bipartition_edge_map[nd.edge.bipartition].head_node.pp
                if nd.pp < population_probability_threshold:
                    new_pop_label = "{}{:03d}".format(
                            self.population_label_prefix,
                            new_pop_idx)
                    new_pop_idx += 1
                    for desc_nd in nd.preorder_iter():
                        if desc_nd.is_leaf():
                            # lineage_population_map[desc_nd.taxon.label] = new_pop_label
                            self.map_sequence_to_pop(
                                    bpp_input_population_label=desc_nd.taxon.label,
                                    bpp_result_population_label=new_pop_label,
                                    lineage_population_map=lineage_population_map)
                            # lineage_population_map[desc_nd.taxon.label] = new_pop_label
                        desc_nd.is_skip_subtree = True
                elif nd.is_leaf():
                    # lineage_population_map[nd.taxon.label] = nd.taxon.label
                    self.map_sequence_to_pop(
                            bpp_input_population_label=nd.taxon.label,
                            bpp_result_population_label=nd.taxon.label,
                            lineage_population_map=lineage_population_map)
        self.logger.info("Posterior probability threshold of {:0.2f}: {} populations".format(population_probability_threshold, len(set(lineage_population_map.values()))))
        return lineage_population_map

    def map_sequence_to_pop(
            self,
            bpp_input_population_label,
            bpp_result_population_label,
            lineage_population_map):
        for seq_label in self.candidate_population_map[bpp_input_population_label]:
            lineage_population_map[seq_label] = bpp_result_population_label

def main():
    """
    Main CLI handler.
    """
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="\n".join(usage_examples),
            )
    parser.add_argument("-i", "--imap",
            action="append",
            nargs="+",
            required=True,
            help="Path to one or more BP&P 'imap' input file(s).")
    parser.add_argument("-r", "--results",
            action="append",
            nargs="+",
            required=True,
            help="Path to one or more BP&P results output file(s).")
    parser.add_argument("-p", "--population-probability-threshold",
            action="append",
            type=float,
            nargs="+",
            default=None,
            help="Probability threshold of population identity.")
    parser.add_argument("--field-delimiter",
            default="\t",
            help="Delimiter for output columns (default=<TAB>)")
    parser.add_argument("-o", "--output-prefix",
            action="store",
            default="coalescent-pops",
            help="Prefix for output files (default=%(default)s).")
    parser.add_argument("-l", "--population-label-prefix",
            action="store",
            default="coalescentpop",
            help="Prefix for delimited (collapsed) population labels (default=%(default)s).")
    parser.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Do not show progress messages.")
    args = parser.parse_args()
    if args.quiet:
        stderr_logging_level = utility.logging.CRITICAL
    else:
        stderr_logging_level = utility.logging.INFO
    bppsum = BppSummarizer(
            population_label_prefix=args.population_label_prefix,
            stderr_logging_level=stderr_logging_level,
            )
    bppsum.read_imap(args.imap)
    bppsum.read_bpp_output(args.results)
    if args.population_probability_threshold is None:
        population_probability_thresholds = (1.0, 0.95, 0.90, 0.75, 0.50)
    else:
        population_probability_thresholds = flatten(args.population_probability_threshold)
    bppsum.calc_populations(
            output_prefix=args.output_prefix,
            population_probability_thresholds=population_probability_thresholds,
            field_delimiter=args.field_delimiter)

if __name__ == '__main__':
    main()

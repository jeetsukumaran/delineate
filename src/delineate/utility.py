#! /usr/bin/env python
# -*- coding: utf-8 -*-

import collections
import re

# CLI Support Functions {{{1

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

def add_source_options(parser, source_options=None):
    if not source_options:
        source_options = parser.add_argument_group("Source options")
    source_options.add_argument("-t", "--tree-file",
            required=True,
            help="Path to tree file.")
    source_options.add_argument("-c", "--config-file",
            help="Path to configuration file.")
    source_options.add_argument("-f", "--format",
            dest="data_format",
            type=str,
            default="nexus",
            choices=["nexus", "newick"],
            help="Input data format (default='%(default)s').")
    return parser

def add_output_options(parser, output_options=None):
    if not output_options:
        output_options = parser.add_argument_group("Output options")
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
    output_options.add_argument( "--field-separator",
            default="\t",
            dest="output_field_separator",
            help="Field separator or delimiter character [default: tab].")
    return parser

# }}}1


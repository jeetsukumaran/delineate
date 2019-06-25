#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import collections
import re
import csv

LINEAGE_ID_FIELDNAME = "lineage"
SPECIES_ID_FIELDNAME = "species"
STATUS_FIELDNAME = "status"

_LOGGING_LEVEL_ENVAR = "DELINEATE_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "DELINEATE_LOGGING_FORMAT"

class RunLogger(object):

    def __init__(self, **kwargs):
        self.name = kwargs.get("name", "RunLog")
        self.is_include_timestamp = kwargs.get("is_include_timestamp", True)
        self.is_include_name = kwargs.get("is_include_name", True)
        self._log = logging.getLogger(self.name)
        self._log.setLevel(logging.DEBUG)
        self.handlers = []
        if kwargs.get("log_to_stderr", True):
            handler1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", logging.INFO))
            handler1.setLevel(stderr_logging_level)
            handler1.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler1)
            self.handlers.append(handler1)
        if kwargs.get("log_to_file", True):
            if "log_stream" in kwargs:
                log_stream = kwargs.get("log_stream")
            else:
                log_stream = open(kwargs.get("log_path", self.name + ".log"), "w")
            handler2 = logging.StreamHandler(log_stream)
            file_logging_level = self.get_logging_level(kwargs.get("file_logging_level", logging.DEBUG))
            handler2.setLevel(file_logging_level)
            handler2.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler2)
            self.handlers.append(handler2)

    def get_logging_level(self, level=None):
        if level in [logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING,
            logging.ERROR, logging.CRITICAL]:
            return level
        elif level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        if level_name == "NOTSET":
            level = logging.NOTSET
        elif level_name == "DEBUG":
            level = logging.DEBUG
        elif level_name == "INFO":
            level = logging.INFO
        elif level_name == "WARNING":
            level = logging.WARNING
        elif level_name == "ERROR":
            level = logging.ERROR
        elif level_name == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
        return level

    def get_default_formatter(self):
        prefix = []
        if self.is_include_name:
            prefix.append(self.name)
        if self.is_include_timestamp:
            prefix.append("%(asctime)s")
        if prefix:
            prefix = "[" + " ".join(prefix) + "] "
        else:
            prefix = ""
        f = logging.Formatter("{}{}".format(prefix, "%(message)s"))
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def debug(self, msg):
        self._log.debug("[DEBUG] {}".format(msg))

    def info(self, msg):
        self._log.info(msg)

    def warning(self, msg):
        self._log.warning(msg)

    def error(self, msg):
        self._log.error(msg)

    def critical(self, msg):
        self._log.critical(msg)

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

def parse_delimited_configuration_file(src,
        delimiter,
        logger):
    src_data = csv.DictReader(src,
        delimiter=delimiter,
        quoting=csv.QUOTE_NONE,
        )
    fieldname_set = set(src_data.fieldnames)
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
    return config_d

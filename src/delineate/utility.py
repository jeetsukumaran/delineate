#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import collections
import re
import csv
import json
import sys

SPECIES_LEAFSET_CONSTRAINTS_KEY = "species_leafset_constraints"
LINEAGE_ID_FIELDNAME = "lineage"
SPECIES_ID_FIELDNAME = "species"
STATUS_FIELDNAME = "status"
CONFIGURATION_REQUIRED_FIELDS = [
        LINEAGE_ID_FIELDNAME,
        SPECIES_ID_FIELDNAME,
        STATUS_FIELDNAME,
    ]

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
    source_options.add_argument("-f", "--tree-format",
            dest="tree_format",
            type=str,
            default="nexus",
            choices=["nexus", "newick"],
            help="Tree file data format (default='%(default)s').")
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

def parse_configuration(args, logger, delimiter=None):
    if args.config_file:
        with open(args.config_file) as src:
            if args.config_file.endswith("json"):
                config = json.load(src)
            else:
                config = parse_delimited_configuration_file(
                        src=src,
                        delimiter=delimiter,
                        logger=logger)
    else:
        config = {}
    return config

def parse_delimited_configuration_file(
        src,
        delimiter=None,
        logger=None):
    if not logger:
        logger = RunLogger(name="delineate",
            is_include_name=True,
            is_include_timestamp=False,
            log_to_stderr=True,
            stderr_logging_level=utility.logging.INFO,
            log_to_file=False,
            file_logging_level=utility.logging.INFO,
            )
    if delimiter is None:
        dialect = csv.Sniffer().sniff(src.read(), delimiters=",\t")
        src.seek(0)
        # delimiter = dialect.delimiter
        src_data = csv.DictReader(src,
                dialect=dialect)
    else:
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
    for required_field in CONFIGURATION_REQUIRED_FIELDS:
        if required_field not in fieldname_set:
            logger.error("ERROR: Required field '{}' not found in configuration source".format(required_field))
            sys.exit(1)
    species_lineage_map = {}
    lineage_species_map = {}
    known = []
    unknown = []
    for entry in src_data:
        if entry[STATUS_FIELDNAME] == "1":
            try:
                species_lineage_map[entry[SPECIES_ID_FIELDNAME]].append(entry[LINEAGE_ID_FIELDNAME])
            except KeyError:
                species_lineage_map[entry[SPECIES_ID_FIELDNAME]] = [entry[LINEAGE_ID_FIELDNAME]]
            known.append(entry[LINEAGE_ID_FIELDNAME])
            lineage_species_map[entry[LINEAGE_ID_FIELDNAME]] = entry[SPECIES_ID_FIELDNAME]
        elif entry[STATUS_FIELDNAME] == "0":
            unknown.append(entry[LINEAGE_ID_FIELDNAME])
            lineage_species_map[entry[LINEAGE_ID_FIELDNAME]] = entry[SPECIES_ID_FIELDNAME]
            pass
        else:
            sys.exit("Unrecognized status: '{}'".format(entry[STATUS_FIELDNAME]))
    # logger.info("{} lineages in total".format(len(known) + len(unknown)))
    # msg = []
    # msg.append("{} species defined in total".format(len(species_lineage_map)))
    # for spidx, sp in enumerate(species_lineage_map):
    #     msg.append("    [{}/{}] '{}'".format(spidx+1, len(species_lineage_map), sp))
    # msg = "\n".join(msg)
    # logger.info(msg)
    # logger.info("{} lineages assigned to {} species".format(len(known), len(species_lineage_map)))
    # logger.info("{} lineages of unknown species affinities".format(len(unknown)))
    species_leafset_constraints = []
    for key in species_lineage_map:
        species_leafset_constraints.append(species_lineage_map[key])
    assert len(species_leafset_constraints) == len(species_lineage_map)
    config_d = {}
    config_d[SPECIES_LEAFSET_CONSTRAINTS_KEY] = species_leafset_constraints
    config_d["configuration_file"] = {}
    config_d["configuration_file"]["lineages"] = known + unknown
    config_d["configuration_file"]["constrained_lineages"] = known
    config_d["configuration_file"]["unconstrained_lineages"] = unknown
    config_d["configuration_file"]["lineage_species_map"] = lineage_species_map
    config_d["configuration_file"]["species_lineage_map"] = species_lineage_map
    return config_d

def report_configuration(
        config_d,
        tree=None,
        logger=None,
        output_file=None,
        output_format="json",
        is_pretty_print=True,
        output_delimiter="\t",
        is_normalize_case=True,
        is_fail_on_extra_tree_lineages=False,
        is_fail_on_extra_configuration_lineages=True,
        ):

    ### VALIDATION and NORMALIZATION

    msg = []
    tree_lineages = None
    config_lineages = None
    if tree is not None:
        tree_lineages = [t.label for t in tree.taxon_namespace]
        msg.append("{} terminal lineages on tree".format(len(tree_lineages)))
    else:
        tree_lineages = None
    if "configuration_file" in config_d and "lineages" in config_d["configuration_file"]:
        config_lineages = list(config_d["configuration_file"]["lineages"])
        msg.append("{} lineages described in configuration file".format(len(config_lineages)))
    else:
        config_lineages = list(lineage_species_map.keys())
    lineage_case_normalization_map = {}
    if tree_lineages is not None and config_lineages is not None:
        tree_lineage_set = set(tree_lineages)
        config_lineage_set = set(config_lineages)
        if is_normalize_case:
            tree_lineages_lower = {}
            for lineage in tree_lineage_set:
                tree_lineages_lower[lineage.lower()] = lineage
                # lineage_case_normalization_map[lineage] = lineage
            for lineage in config_lineage_set:
                if lineage.lower() in tree_lineages_lower:
                    lineage_case_normalization_map[lineage] = tree_lineages_lower[lineage.lower()]
                else:
                    lineage_case_normalization_map[lineage] = lineage
            config_lineage_set_normalized = set([lineage_case_normalization_map[lineage] for lineage in config_lineage_set])
            config_lineage_set = config_lineage_set_normalized
        s1 = tree_lineage_set - config_lineage_set
        if s1:
            msg.append("{} lineages in tree not described in configuration file: {}".format(len(s1), ", ".join(s1)))
            s1_error_msg = ["{}: {} lineages found on tree but not in configuration file:".format(
                "ERROR" if is_fail_on_extra_tree_lineages else "NOTE",
                len(s1))]
            for lidx, label in enumerate(sorted(s1, key=lambda x: x.lower())):
                s1_error_msg.append("    [{: 3d}/{:<3d}] {}".format(lidx+1, len(s1), label))
            s1_error_msg = "\n".join(s1_error_msg)
        else:
            s1_error_msg = ""
        s2 = config_lineage_set - tree_lineage_set
        if s2:
            msg.append("{} lineages in configuration file not found on tree: {}".format(len(s2), ", ".join(s2)))
            s2_error_msg = ["{}: {} lineages found in configuration but not found on tree:".format(
                "ERROR" if is_fail_on_extra_tree_lineages else "WARNING",
                len(s2))]
            for lidx, label in enumerate(sorted(s2, key=lambda x: x.lower())):
                s2_error_msg.append("    [{: 3d}/{:<3d}] {}".format(lidx+1, len(s2), label))
            s2_error_msg = "\n".join(s2_error_msg)
        else:
            s2_error_msg = ""
        if s1 or s2:
            if s1 and is_fail_on_extra_tree_lineages:
                logger.error(s1_error_msg)
            elif s1:
                logger.warning(s1_error_msg)
            if s2 and is_fail_on_extra_configuration_lineages:
                logger.error(s2_error_msg)
            elif s2:
                logger.warning(s2_error_msg)
            if s1 and is_fail_on_extra_tree_lineages:
                logger.error("Exiting due to lineage identity mismatch error")
                sys.exit(1)
            if s2 and is_fail_on_extra_configuration_lineages:
                logger.error("Exiting due to lineage identity mismatch error")
                sys.exit(1)
        all_lineages = tree_lineages
    elif tree_lineages is not None:
        all_lineages = tree_lineages
    elif config_lineages is not None:
        all_lineages = config_lineages

    ### MAPPING

    species_lineage_map = {}
    lineage_species_map = {}
    original_lineage_species_map = {}
    if "configuration_file" in config_d:
        # additional info from configurator front-end
        for spp in config_d["configuration_file"]["species_lineage_map"]:
            species_lineage_map[spp] = []
            for lineage_name in config_d["configuration_file"]["species_lineage_map"][spp]:
                species_lineage_map[spp].append(lineage_case_normalization_map.get(lineage_name, lineage_name))
        lineage_species_map = {}
        original_lineage_species_map.update(config_d["configuration_file"]["lineage_species_map"])
        for lineage_name in config_d["configuration_file"]["lineage_species_map"]:
            normalized_lineage_name = lineage_case_normalization_map.get(lineage_name, lineage_name)
            if lineage_name in config_d["configuration_file"]["constrained_lineages"]:
                lineage_species_map[normalized_lineage_name] = config_d["configuration_file"]["lineage_species_map"][lineage_name]
    else:
        species_leafsets = config_d[SPECIES_LEAFSET_CONSTRAINTS_KEY]
        if isinstance(species_leafsets, list) or isinstance(species_leafsets, tuple):
            for spi, sp in enumerate(species_leafsets):
                lineages = []
                for lineage in sp:
                    normalized_lineage_name = lineage_case_normalization_map.get(lineage, lineage)
                    lineage_species_map[normalized_lineage_name] = sp_label
                    lineages.append(normalized_lineage_name)
                sp_label = "Sp{:03d}".format(spi+1)
                species_lineage_map[sp_label] = lineages
        else:
            raise NotImplementedError
            species_lineage_map = {}
            for sp_label in species_leafsets:
                for lineage in species_leafsets[sp_label]:
                    lineages = []
                    normalized_lineage_name = lineage_case_normalization_map.get(lineage, lineage)
                    lineage_species_map[normalized_lineage_name] = sp_label
                    lineages.append(normalized_lineage_name)
                species_lineage_map[sp_label] = lineages

    msg.append("{} lineages with known species identities, assigned to {} species".format(
        len(lineage_species_map),
        len(species_lineage_map)))
    msg.append("{} lineages with species identities to be inferred".format(
        len(all_lineages) - len(lineage_species_map)))
    if logger:
        pmsg = "\n".join(["  -  {}".format(m) for m in msg])
        logger.info("Analysis configuration:\n{}".format(pmsg))
    if output_file is None:
        pass
    elif output_format == "json":
        kwargs = {}
        if is_pretty_print:
            kwargs["sort_keys"] = True
            kwargs["indent"] = 4
        json.dump(config_d, output_file, **kwargs)
    elif output_format == "delimited":
        output_file.write("{}\n".format(output_delimiter.join(CONFIGURATION_REQUIRED_FIELDS)))
        for lineage in all_lineages:
            parts = []
            parts.append(lineage)
            if lineage in lineage_species_map:
                parts.append(lineage_species_map[lineage])
                parts.append("1")
            else:
                parts.append(original_lineage_species_map.get(lineage, "?"))
                parts.append("0")
            output_file.write("{}\n".format(output_delimiter.join(parts)))
    else:
        raise ValueError(output_format)
    return msg


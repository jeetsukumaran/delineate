#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import collections
import re
import csv
import json
import os
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
    config_d["configuration_file"]["constrained_lineage_species_map"] = {lineage_name:lineage_species_map[lineage_name] for lineage_name in known}
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
        msg.append("{} terminal lineages on population tree".format(len(tree_lineages)))
    else:
        tree_lineages = None
    if "configuration_file" in config_d and "lineages" in config_d["configuration_file"]:
        config_lineages = list(config_d["configuration_file"]["lineages"])
        msg.append("{} lineages described in configuration file".format(len(config_lineages)))
    else:
        if tree_lineages is not None:
            config_lineages = list(tree_lineages)
        else:
            raise ValueError("Tree file required for checking JSON file configuration")
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
            # msg.append("NOTE: {} lineages in tree not described in configuration file: {}".format(len(s1), ", ".join(s1)))
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
            # msg.append("WARNING: {} lineages in configuration file not found on tree: {}".format(len(s2), ", ".join(s2)))
            s2_error_msg = ["{}: {} lineages described in configuration file but not found on tree:".format(
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
                msg.append(s1_error_msg)
            if s2 and is_fail_on_extra_configuration_lineages:
                logger.error(s2_error_msg)
            elif s2:
                msg.append(s2_error_msg)
            if s1 and is_fail_on_extra_tree_lineages:
                logger.error("Exiting due to lineage identity mismatch error")
                sys.exit(1)
            if s2 and is_fail_on_extra_configuration_lineages:
                logger.error("Exiting due to lineage identity mismatch error")
                sys.exit(1)
        all_lineages = sorted(set(tree_lineages + config_lineages))
    elif tree_lineages is not None:
        all_lineages = sorted(tree_lineages)
    elif config_lineages is not None:
        all_lineages = sorted(config_lineages)

    ### MAPPING

    species_lineage_map = {}
    constrained_lineage_species_map = {}
    full_lineage_species_map = {}
    seen_lineages = set()
    if "configuration_file" in config_d:
        # additional info from configurator front-end
        for spp in config_d["configuration_file"]["species_lineage_map"]:
            species_lineage_map[spp] = []
            for lineage_name in config_d["configuration_file"]["species_lineage_map"][spp]:
                normalized_lineage_name = lineage_case_normalization_map.get(lineage_name, lineage_name)
                if normalized_lineage_name in seen_lineages:
                    logger.error("ERROR: Duplicate lineage species assignment: '{}'".format(normalized_lineage_name))
                    sys.exit(1)
                seen_lineages.add(normalized_lineage_name)
                species_lineage_map[spp].append(normalized_lineage_name)
        constrained_lineage_species_map = {}
        for lineage_name in config_d["configuration_file"]["lineage_species_map"]:
            normalized_lineage_name = lineage_case_normalization_map.get(lineage_name, lineage_name)
            full_lineage_species_map[normalized_lineage_name] = config_d["configuration_file"]["lineage_species_map"][lineage_name]
            if lineage_name in config_d["configuration_file"]["constrained_lineages"]:
                constrained_lineage_species_map[normalized_lineage_name] = full_lineage_species_map[normalized_lineage_name]
    elif SPECIES_LEAFSET_CONSTRAINTS_KEY in config_d:
        species_leafsets = config_d[SPECIES_LEAFSET_CONSTRAINTS_KEY]
        if isinstance(species_leafsets, list) or isinstance(species_leafsets, tuple):
            for spi, sp in enumerate(species_leafsets):
                lineages = []
                sp_label = compose_constrained_species_label(spi)
                for lineage in sp:
                    normalized_lineage_name = lineage_case_normalization_map.get(lineage, lineage)
                    constrained_lineage_species_map[normalized_lineage_name] = sp_label
                    lineages.append(normalized_lineage_name)
                species_lineage_map[sp_label] = lineages
        else:
            raise NotImplementedError
            species_lineage_map = {}
            for sp_label in species_leafsets:
                for lineage in species_leafsets[sp_label]:
                    lineages = []
                    normalized_lineage_name = lineage_case_normalization_map.get(lineage, lineage)
                    constrained_lineage_species_map[normalized_lineage_name] = sp_label
                    lineages.append(normalized_lineage_name)
                species_lineage_map[sp_label] = lineages

    if species_lineage_map:
        spp_list = []
        species_names = species_lineage_map.keys()
        max_spp_name_length = max(len(sp) for sp in species_names)
        sp_name_template = "{{:{}}}".format(max_spp_name_length + 2)
        spp_list.append("{} species defined in configuration file:".format(
            len(species_lineage_map)))
        for sidx, spp in enumerate(sorted(species_names)):
            if len(species_lineage_map[spp]) > 1:
                descriptor = "lineages"
            else:
                descriptor = "lineage"
            spp_list.append("    [{: 3d}/{:<3d}] SPECIES: {} ({} {})".format(
                    sidx+1,
                    len(species_lineage_map),
                    sp_name_template.format("'"+spp+"'"),
                    len(species_lineage_map[spp]),
                    descriptor))
        msg.append("\n".join(spp_list))
    else:
        msg.append("0 species defined in configuration file")

    constrained_lineage_list =[]
    constrained_lineage_list.append("{} lineages with known species affinities:".format(
        len(constrained_lineage_species_map),
        # len(species_lineage_map)),
        ))
    max_lineage_name_length = max(len(ln) for ln in all_lineages)
    lineage_name_template = "{{:{}}}".format(max_lineage_name_length + 2)
    lidx = 0
    for lineage in all_lineages:
        if lineage in constrained_lineage_species_map:
            lidx += 1
            constrained_lineage_list.append("    [{: 3d}/{:<3d}] LINEAGE {} (SPECIES: '{}')".format(
                lidx,
                len(all_lineages) - len(constrained_lineage_species_map),
                lineage_name_template.format("'"+lineage+"'"),
                constrained_lineage_species_map[lineage]
                ))
    msg.append("\n".join(constrained_lineage_list))

    unconstrained_lineage_list =[]
    unconstrained_lineage_list.append("{} lineages of unknown species affinities:".format(
        len(all_lineages) - len(constrained_lineage_species_map)))
    lidx = 0
    for lineage in all_lineages:
        if lineage not in constrained_lineage_species_map:
            lidx += 1
            unconstrained_lineage_list.append("    [{: 3d}/{:<3d}] LINEAGE '{}'".format(
                lidx,
                len(all_lineages) - len(constrained_lineage_species_map),
                lineage))
    msg.append("\n".join(unconstrained_lineage_list))
    if logger:
        pmsg = "\n".join(["  - {}".format(m) for m in msg])
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
        for lidx, lineage in enumerate(all_lineages):
            parts = []
            parts.append(lineage)
            if lineage in constrained_lineage_species_map:
                parts.append(constrained_lineage_species_map[lineage])
                parts.append("1")
            else:
                parts.append(full_lineage_species_map.get(lineage, "?"))
                parts.append("0")
            output_file.write("{}\n".format(output_delimiter.join(parts)))
    else:
        raise ValueError(output_format)
    for lineage_name in lineage_case_normalization_map:
        if lineage_name not in constrained_lineage_species_map:
            normalized_lineage_name = lineage_case_normalization_map.get(lineage_name, None)
            if normalized_lineage_name in constrained_lineage_species_map:
                constrained_lineage_species_map[lineage_name] = constrained_lineage_species_map[normalized_lineage_name]
    return {
            "constrained_lineage_species_map": constrained_lineage_species_map,
            "full_lineage_species_map": full_lineage_species_map,
            "species_lineage_map": species_lineage_map,
            "lineage_case_normalization_map": lineage_case_normalization_map,
            "messages": msg,
        }

def compose_constrained_species_label(idx):
    return "ConstrainedSp{:03d}".format(idx+1)

def compose_output_prefix(input_filepath, default):
    if not input_filepath:
        return default
    return os.path.splitext(os.path.expanduser(os.path.expandvars(input_filepath)))[0]

def compose_lineage_species_name_map(
        postanalysis_species_leafset_labels,
        preanalysis_lineage_species_name_map,
        is_validate_species_group_consistency=True,
        ):
    lnsp_map = {}
    new_sp_idx = 0
    existing_sp_names = set(preanalysis_lineage_species_name_map.values())
    processed_sp_names = set()
    for leafset in postanalysis_species_leafset_labels:
        sp_id = None
        for lineage in leafset:
            try:
                assigned_sp_id = preanalysis_lineage_species_name_map[lineage]
                if (is_validate_species_group_consistency and
                    (sp_id is not None and assigned_sp_id != sp_id)):
                    msg = []
                    msg.append("Conflicting species identity assignment for lineages grouped into same species '{}' vs '{}':".format(
                        sp_id,
                        assigned_sp_id))
                    max_lineage_name_length = max(len(ln) for ln in all_lineages)
                    lineage_name_template = "{{:{}}}".format(max_lineage_name_length + 2)
                    for lineage in leafset:
                        msg.append("  - LINEAGE {} => SPECIES {}".format(
                            lineage_name_template.format("'{}'".format(lineage)),
                            assigned_sp_id))
                    msg = "\n".join(msg)
                    raise ValueError(msg)
                sp_id = assigned_sp_id
                # break
            except KeyError:
                pass
        if sp_id is None:
            while True:
                new_sp_idx += 1
                sp_id = "DelineatedSp{:03d}".format(new_sp_idx)
                if sp_id not in existing_sp_names:
                    existing_sp_names.add(sp_id)
                    break
        assert sp_id is not None
        assert sp_id not in processed_sp_names
        processed_sp_names.add(sp_id)
        for lineage in leafset:
            lnsp_map[lineage] = sp_id
    return lnsp_map


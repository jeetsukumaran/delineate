#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import collections
import re
import csv
import json
import os
import sys

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
        if kwargs.get("is_log_to_stderr", True):
            handler1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", logging.INFO))
            handler1.setLevel(stderr_logging_level)
            handler1.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler1)
            self.handlers.append(handler1)
        if kwargs.get("is_log_to_file", True):
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

def compose_table(**kwargs):
    return "\n".join(compose_table_rows(**kwargs))

def compose_table_rows(
        columns,
        prefixes,
        quoted,
        is_indexed=True,
        indent="",
        ):
    s = []
    field_templates = []
    for column, is_quoted in zip(columns, quoted):
        ml = 0
        for value in column:
            ml = max(ml, len(value))
        if is_quoted:
            ml += 2
        field_templates.append("{{:{}}}".format(ml))
        # field_templates.append("{:10}")
    for row_idx in range(len(columns[0])):
        row = []
        if indent:
            row.append(indent)
        if is_indexed:
            row.append("[{: 3d}/{:<3d}]".format(row_idx+1, len(columns[0])))
        for column, prefix, field_template, is_quoted in zip(columns, prefixes, field_templates, quoted):
            if prefix:
                row.append(prefix)
            if is_quoted:
                field_text = "'{}'".format(column[row_idx])
            else:
                field_text = column[row_idx]
            row.append(" ")
            row.append(field_template.format(field_text))
        s.append("".join(row))
    return s


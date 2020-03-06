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

class ColorMap(object):
    """
    Colors that contrast well over various color perception regimes.
    """

    contrast_pairs = [
            ["#ffc20a", "#0c7bdc"],
            ["#1aff1a", "#4b0092"],
            ["#994f00", "#006cd1"],
            ["#fefe62", "#d35fb7"],
            ["#e1be6a", "#40b0a6"],
            ["#005ab5", "#dc3220"],
            ["#e66100", "#5d3a9b"],
            ["#1a85ff", "#d41159"],
    ]

    def __init__(self):
        self.label_color_map = {}
        self.color_label_map = {}
        self.colors = [
            "#000000",
            "#e69f00",
            "#56b4e9",
            "#009e73",
            "#f0e442",
            "#0072b2",
            "#d55e00",
            "#cc79a7"
        ]
        self.available_colors = list(self.colors)

    def __call__(self, label):
        try:
            return self.label_color_map[label]
        except KeyError:
            new_color = self.available_colors.pop()
            self.label_color_map[label] = new_color
            self.color_label_map[new_color] = label
        return new_color

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

def error_exit(msg, logger):
    logger.error("ERROR: {}".format(msg))
    logger.critical("Terminating due to error")
    sys.exit(1)

def compose_output_prefix(input_filepath, default):
    if not input_filepath:
        return default
    return os.path.basename(os.path.expanduser(os.path.expandvars(input_filepath)))

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
        if ml > 0:
            field_templates.append("{{:{}}}".format(ml))
        else:
            field_templates.append("")
        # field_templates.append("{:10}")
    for row_idx in range(len(columns[0])):
        row = []
        if indent:
            row.append(indent)
        if is_indexed:
            row.append("[{: 3d}/{:<3d}] ".format(row_idx+1, len(columns[0])))
        for col_idx, (column, prefix, field_template, is_quoted) in enumerate(zip(columns, prefixes, field_templates, quoted)):
            if col_idx > 0:
                row.append("    ")
            if prefix:
                row.append(prefix)
            if is_quoted:
                field_text = "'{}'".format(column[row_idx])
            else:
                field_text = column[row_idx]
            row.append(field_template.format(field_text))
        s.append("".join(row))
    return s


#! /usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import csv
import json
import sys
from delineate import model
from delineate import utility
from delineate import summarize
from dendropy.utility.container import OrderedCaselessDict

SPECIES_LEAFSET_CONSTRAINTS_KEY = "species_leafset_constraints"
LINEAGE_ID_FIELDNAME = "lineage"
SPECIES_ID_FIELDNAME = "species"
STATUS_FIELDNAME = "status"
CONFIGURATION_REQUIRED_FIELDS = [
        LINEAGE_ID_FIELDNAME,
        SPECIES_ID_FIELDNAME,
        STATUS_FIELDNAME,
    ]

def get_controller(
        args,
        name,
        logger=None,
        logger_kwargs=None):
    controller = Controller(
            name=name,
            is_case_sensitive=getattr(args, "is_case_sensitive", False),
            is_fail_on_extra_tree_lineages=getattr(args, "is_fail_on_extra_tree_lineages", True),
            is_fail_on_extra_configuration_lineages=getattr(args, "is_fail_on_extra_configuration_lineages", True),
            logger=logger,
            logger_kwargs=logger_kwargs)
    controller.read_tree(
            tree_filepath=args.tree_file,
            schema=args.tree_format,
            underflow_protection=getattr(args, "underflow_protection", False))
    controller.parse_configuration_file(
            config_filepath=args.config_file,
            delimiter=None)
    controller.register_names()
    return controller

class Registry(object):

    def __init__(self, **kwargs):
        self.tree_lineage_names = None
        self.config_lineage_names = None
        self.is_case_sensitive = kwargs.pop("is_case_sensitive", False)
        self.is_fail_on_extra_tree_lineages = kwargs.pop("is_fail_on_extra_tree_lineages", True)
        self.is_fail_on_extra_configuration_lineages = kwargs.pop("is_fail_on_extra_configuration_lineages", True)
        self.logger = kwargs.pop("logger")
        self.original_to_normalized_lineage_name_map = {}
        self.config_name_normalization_report = {}
        self.preanalysis_constrained_species_lineages_map = {}
        if self.is_case_sensitive:
            self.normalized_tree_lineage_names = {}
            self.normalized_config_lineage_names = {}
            self.species_names = {}
            self.preanalysis_constrained_lineage_species_map = {}
        else:
            self.normalized_tree_lineage_names = OrderedCaselessDict()
            self.normalized_config_lineage_names = OrderedCaselessDict()
            self.normalized_species_names = OrderedCaselessDict()
            self.preanalysis_constrained_lineage_species_map = OrderedCaselessDict()
        self.extra_tree_lineage_names = []
        self.extra_configuration_lineages = []

    def normalize_lineage_names(self):
        tree_lineage_set = set(self.tree_lineage_names)
        # tree lineages give the canonical orthography
        for lineage in self.tree_lineage_names:
            self.normalized_tree_lineage_names[lineage] = lineage
            self.original_to_normalized_lineage_name_map[lineage] = lineage
        normalized_configuration_lineages = {}
        extra_configuration_lineages = set()
        for lineage in self.config_lineage_names:
            self.normalized_config_lineage_names[lineage] = lineage
            try:
                normalized_name = self.normalized_tree_lineage_names[lineage]
                self.original_to_normalized_lineage_name_map[lineage] = normalized_name
                if normalized_name != lineage:
                    self.config_name_normalization_report[lineage] = "(NORMALIZED TO: '{}')".format(normalized_name)
                    normalized_configuration_lineages[lineage] = normalized_name
                else:
                    self.config_name_normalization_report[lineage] = ""
            except KeyError as e:
                # This is a serious error: it means that the configuration file
                # has a taxon that is not on the tree. But we handle this issue
                # later so a full report can be shown
                self.config_name_normalization_report[lineage] = "(NOT FOUND ON TREE)"
                extra_configuration_lineages.add(lineage)
        self.normalization_report(
                normalized_configuration_lineages=normalized_configuration_lineages,
                extra_configuration_lineages=extra_configuration_lineages)

    def read_configuration_table_species(self,
            conf_lineage_species_map,
            conf_constrained_lineages):
        if self.is_case_sensitive:
            nccl = {}
        else:
            nccl = OrderedCaselessDict()
        for ln in conf_constrained_lineages:
            nccl[ln] = True
        for lineage_name in conf_lineage_species_map:
            if lineage_name not in nccl:
                continue
            species_name = conf_lineage_species_map[lineage_name]
            if species_name not in self.normalized_species_names:
                self.normalized_species_names[species_name] = species_name
            else:
                species_name = self.normalized_species_names[species_name]
            try:
                normalized_lineage_name = self.original_to_normalized_lineage_name_map[lineage_name]
            except KeyError:
                utility.error_exit(
                        msg="Lineage '{}' not defined (missing on tree?)".format(lineage_name),
                        logger=self.logger)
            if normalized_lineage_name in self.preanalysis_constrained_lineage_species_map:
                utility.error_exit(
                        msg="Duplicate lineage species assignment: '{}'".format(normalized_lineage_name),
                        logger=self.logger)
            self.preanalysis_constrained_lineage_species_map[normalized_lineage_name] = species_name
            try:
                self.preanalysis_constrained_species_lineages_map[species_name].add(normalized_lineage_name)
            except KeyError:
                self.preanalysis_constrained_species_lineages_map[species_name] = set([normalized_lineage_name])
        self.preanalysis_constrained_species_report()

    def compile_configuration_species_groupings(self, species_leafset_constraints):
        for spi, sp in enumerate(species_leafset_constraints):
            lineages = []
            species_name = "ConstrainedSp{:03d}".format(spi+1)
            self.normalized_species_names[species_name] = species_name
            for lineage_name in sp:
                try:
                    normalized_lineage_name = self.original_to_normalized_lineage_name_map[lineage_name]
                except KeyError:
                    utility.error_exit(
                            msg="Lineage '{}' not defined (missing on tree?)".format(lineage_name),
                            logger=self.logger)
                self.preanalysis_constrained_lineage_species_map[normalized_lineage_name] = species_name
                try:
                    self.preanalysis_constrained_species_lineages_map[species_name].add(normalized_lineage_name)
                except KeyError:
                    self.preanalysis_constrained_species_lineages_map[species_name] = set([normalized_lineage_name])
        self.preanalysis_constrained_species_report()

    def preanalysis_constrained_species_report(self):
        species_names = sorted(self.preanalysis_constrained_species_lineages_map.keys())
        num_lineages = ["({} lineages)".format(len(self.preanalysis_constrained_species_lineages_map[n])) for n in species_names]
        stbl = utility.compose_table(
                columns=[
                    species_names,
                    num_lineages,
                    ],
                prefixes=["", ""],
                quoted=[True, False],
                is_indexed=True,
                indent="    ")
        self.logger.info("{} species defined in configuration constraints, with {} lineages assigned:\n{}".format(
                len(species_names),
                len(self.preanalysis_constrained_lineage_species_map),
                stbl,
                ))
        constrained_lineages = sorted(self.preanalysis_constrained_lineage_species_map.keys(), key=lambda n: (self.preanalysis_constrained_lineage_species_map[n], n))
        species_assignments = ["(SPECIES: '{}')".format(self.preanalysis_constrained_lineage_species_map[n]) for n in constrained_lineages]
        lntbl = utility.compose_table(
                columns=[
                    constrained_lineages,
                    species_assignments,
                    ],
                prefixes=["", ""],
                quoted=[True, False],
                is_indexed=True,
                indent="    ")
        self.logger.info("{} out of {} lineages assigned by constraints to {} species:\n{}".format(
                len(constrained_lineages),
                len(self.tree_lineage_names),
                len(species_names),
                lntbl,
                ))
        unconstrained_lineages = sorted(n for n in self.tree_lineage_names if n not in self.preanalysis_constrained_lineage_species_map)
        lntbl = utility.compose_table(
                columns=[
                    unconstrained_lineages,
                    ],
                prefixes=[""],
                quoted=[True],
                is_indexed=True,
                indent="    ")
        self.logger.info("{} out of {} lineages not constrained by species assignments:\n{}".format(
                len(unconstrained_lineages),
                len(self.tree_lineage_names),
                lntbl,
                ))
        assert len(unconstrained_lineages) + len(constrained_lineages) == len(self.tree_lineage_names)

    def normalization_report(self,
            normalized_configuration_lineages,
            extra_configuration_lineages):
        treetbl = utility.compose_table(
                columns=[self.tree_lineage_names,
                    ["(NOT FOUND IN CONFIGURATION)" if lineage not in self.normalized_config_lineage_names else "" for lineage in self.tree_lineage_names],
                    ],
                prefixes=["", ""],
                quoted=[True, False],
                is_indexed=True,
                indent="    ")
        self.logger.info("{} lineages found on population tree:\n{}".format(
            len(self.tree_lineage_names),
            treetbl,
            ))
        if extra_configuration_lineages:
            cfntbl = utility.compose_table(
                    columns=[self.config_lineage_names,
                        [self.config_name_normalization_report[n] for n in self.config_lineage_names]
                        ],
                    prefixes=["", ""],
                    quoted=[True, False],
                    is_indexed=True,
                    indent="    ")
            self.logger.info("{} lineages found in configuration file:\n{}".format(
                len(self.config_lineage_names),
                cfntbl,
                ))
        elif normalized_configuration_lineages:
            n1 = list(normalized_configuration_lineages.keys())
            n2 = [normalized_configuration_lineages[k] for k in n1]
            cfntbl = utility.compose_table(
                    columns=[n1, n2, ],
                    prefixes=["", "NORMALIZED TO: "],
                    quoted=[True, True],
                    is_indexed=True,
                    indent="    ")
            self.logger.info("{} lineages found in configuration file, with the following normalized for concordance with tree lineages:\n{}".format(
                len(self.config_lineage_names),
                cfntbl,
                ))
        else:
            self.logger.info("{} lineages found in configuration file fully concordant with tree lineages".format(
                len(self.config_lineage_names),
                ))

    def validate_lineage_names(self):
        for lineage in self.config_lineage_names:
            if lineage not in self.normalized_tree_lineage_names:
                self.extra_configuration_lineages.append(lineage)
        for lineage in self.tree_lineage_names:
            if lineage not in self.normalized_config_lineage_names:
                self.extra_tree_lineage_names.append(lineage)
        if self.extra_tree_lineage_names:
            s1_error_msg = ["{}: {} lineages found on tree but not in configuration data:".format(
                "ERROR" if self.is_fail_on_extra_tree_lineages else "WARNING",
                len(self.extra_tree_lineage_names))]
            s1_error_msg.append(self.compose_name_list(self.extra_tree_lineage_names))
            s1_error_msg = "\n".join(s1_error_msg)
        else:
            s1_error_msg = ""

        if self.extra_configuration_lineages:
            s2_error_msg = ["{}: {} lineages found in configuration data but not on tree:".format(
                "ERROR" if self.is_fail_on_extra_configuration_lineages else "WARNING",
                len(self.extra_configuration_lineages))]
            s2_error_msg.append(self.compose_name_list(self.extra_configuration_lineages))
            s2_error_msg = "\n".join(s2_error_msg)
        else:
            s2_error_msg = ""

        is_fail = []
        if self.extra_tree_lineage_names and self.is_fail_on_extra_tree_lineages:
            self.logger.error(s1_error_msg)
            is_fail.append("1")
        elif s1_error_msg:
            self.logger.warning(s1_error_msg)
        if self.extra_configuration_lineages and self.is_fail_on_extra_configuration_lineages:
            self.logger.error(s2_error_msg)
            is_fail.append("2")
        elif s2_error_msg:
            self.logger.warning(s2_error_msg)
        if is_fail:
            utility.error_exit(
                msg="Lineage identity errors found ({})".format(", ".join(is_fail)),
                logger=self.logger)

    def compose_name_list(self, names):
        s = utility.compose_table(
                columns=[names],
                prefixes=[""],
                quoted=[True],
                is_indexed=True,
                indent="    ")
        return s

    def compose_report(self):
        msg = []
        msg.append("{} terminal lineages on population tree".format(len(self.tree_lineage_names)))
        msg.append("{} lineages described in configuration file".format(len(self.config_lineage_names)))

class Controller(object):

    def __init__(self,
            name,
            is_case_sensitive=True,
            is_fail_on_extra_tree_lineages=False,
            is_fail_on_extra_configuration_lineages=True,
            logger=None,
            logger_kwargs=None):
        self.name = name
        self.is_case_sensitive = is_case_sensitive
        self.is_fail_on_extra_tree_lineages = is_fail_on_extra_tree_lineages
        self.is_fail_on_extra_configuration_lineages = is_fail_on_extra_configuration_lineages
        self._logger = logger
        if logger_kwargs:
            self._logger_kwargs = dict(logger_kwargs)
        else:
            self._logger_kwargs = {}
        self._speciation_completion_rate = None
        self._species_leafset_constraint_labels = None
        self._constrained_lineage_leaf_labels = None
        self.registry = Registry(
                is_case_sensitive=self.is_case_sensitive,
                is_fail_on_extra_tree_lineages=self.is_fail_on_extra_tree_lineages,
                is_fail_on_extra_configuration_lineages=self.is_fail_on_extra_configuration_lineages,
                logger=self.logger,
                )
        self.config_d = {}
        self.tree = None

    @property
    def logger(self):
        if self._logger is None:
            self._logger = utility.RunLogger(
                name=self._logger_kwargs.pop("name", self.name),
                is_include_name=self._logger_kwargs.pop("is_include_name", True),
                is_include_timestamp=self._logger_kwargs.pop("is_include_timestamp", True),
                is_log_to_stderr=self._logger_kwargs.pop("is_log_to_stderr", True),
                stderr_logging_level=utility.logging.INFO,
                is_log_to_file=self._logger_kwargs.pop("is_log_file", False),
                file_logging_level=utility.logging.INFO,
                **self._logger_kwargs,
                )
        return self._logger
    @logger.setter
    def logger(self, value):
        self._logger = value
    @logger.deleter
    def logger(self):
        del self._logger

    @property
    def speciation_completion_rate(self):
        if self._speciation_completion_rate is None:
            try:
                self._speciation_completion_rate = self.config_d["speciation_completion_rate"]
            except KeyError as e:
                return None
        return self._speciation_completion_rate
    @speciation_completion_rate.setter
    def speciation_completion_rate(self, value):
        self._speciation_completion_rate = value
    @speciation_completion_rate.deleter
    def speciation_completion_rate(self):
        del self._speciation_completion_rate

    @property
    def species_leafset_constraint_labels(self):
        if self._species_leafset_constraint_labels is None:
            try:
                self._species_leafset_constraint_labels = self.config_d[SPECIES_LEAFSET_CONSTRAINTS_KEY]
            except KeyError as e:
                pass
        return self._species_leafset_constraint_labels
    @species_leafset_constraint_labels.setter
    def species_leafset_constraint_labels(self, value):
        self._species_leafset_constraint_labels = value
    @species_leafset_constraint_labels.deleter
    def species_leafset_constraint_labels(self):
        del self._species_leafset_constraint_labels

    @property
    def has_species_constraints(self):
        return SPECIES_LEAFSET_CONSTRAINTS_KEY in self.config_d

    @property
    def constrained_lineage_leaf_labels(self):
        if self._constrained_lineage_leaf_labels is None:
            self._constrained_lineage_leaf_labels = list(itertools.chain(*self.species_leafset_constraint_labels))
        return self._constrained_lineage_leaf_labels

    def read_tree(self,
            tree_filepath,
            schema,
            underflow_protection=True):
        self.tree = model.LineageTree.get(
                path=tree_filepath,
                schema=schema,
                )
        self.tree.is_use_decimal_value_type = underflow_protection
        return self.tree

    def parse_configuration_file(self,
            config_filepath,
            delimiter=None,
            ):
        if config_filepath:
            with open(config_filepath) as src:
                if config_filepath.endswith("json"):
                    self.config_d = json.load(src)
                else:
                    self.config_d = self.parse_configuration_table_file(
                            src=src,
                            delimiter=delimiter)
        else:
            self.config_d = {}
        return self.config_d

    def parse_configuration_table_file(
            self,
            src,
            delimiter=None,
            ):
        self.config_d = {}
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
        for required_field in CONFIGURATION_REQUIRED_FIELDS:
            if required_field not in fieldname_set:
                utility.error_exit(
                        msg="Required field '{}' not found in configuration source".format(required_field),
                        logger=logger)
        species_constrained_lineage_map = {}
        lineage_species_map = {}
        known = []
        unknown = []
        for entry in src_data:
            if entry[STATUS_FIELDNAME] == "1":
                try:
                    species_constrained_lineage_map[entry[SPECIES_ID_FIELDNAME]].append(entry[LINEAGE_ID_FIELDNAME])
                except KeyError:
                    species_constrained_lineage_map[entry[SPECIES_ID_FIELDNAME]] = [entry[LINEAGE_ID_FIELDNAME]]
                known.append(entry[LINEAGE_ID_FIELDNAME])
                lineage_species_map[entry[LINEAGE_ID_FIELDNAME]] = entry[SPECIES_ID_FIELDNAME]
            elif entry[STATUS_FIELDNAME] == "0":
                unknown.append(entry[LINEAGE_ID_FIELDNAME])
                lineage_species_map[entry[LINEAGE_ID_FIELDNAME]] = entry[SPECIES_ID_FIELDNAME]
                pass
            else:
                utility.error_exit(
                        msg="Unrecognized status: '{}'".format(entry[STATUS_FIELDNAME]),
                        logger=self.logger)
        species_leafset_constraints = []
        for key in species_constrained_lineage_map:
            species_leafset_constraints.append(species_constrained_lineage_map[key])
        assert len(species_leafset_constraints) == len(species_constrained_lineage_map)
        self.config_d = {}
        self.config_d[SPECIES_LEAFSET_CONSTRAINTS_KEY] = species_leafset_constraints
        self.config_d["configuration_table"] = {}
        self.config_d["configuration_table"]["lineages"] = known + unknown
        self.config_d["configuration_table"]["constrained_lineages"] = known
        self.config_d["configuration_table"]["unconstrained_lineages"] = unknown
        self.config_d["configuration_table"]["lineage_species_map"] = lineage_species_map
        self.config_d["configuration_table"]["constrained_lineage_species_map"] = {lineage_name:lineage_species_map[lineage_name] for lineage_name in known}
        self.config_d["configuration_table"]["species_constrained_lineage_map"] = species_constrained_lineage_map
        return self.config_d

    def register_names(self):
        self.register_lineage_names()
        self.register_preanalysis_species_names()

    def register_lineage_names(self):
        if self.tree is None:
            raise ValueError("'tree' not set")
        self.registry.tree_lineage_names = [t.label for t in self.tree.taxon_namespace]
        if "configuration_table" in self.config_d and "lineages" in self.config_d["configuration_table"]:
            self.registry.config_lineage_names = list(self.config_d["configuration_table"]["lineages"])
        else:
            self.registry.config_lineage_names = list(self.registry.tree_lineage_names)
            # self.registry.config_lineage_names = []
        if not self.registry.config_lineage_names:
            return
        self.registry.normalize_lineage_names()
        self.registry.validate_lineage_names()

    def register_preanalysis_species_names(self):
        species_constrained_lineage_map = {}
        constrained_lineage_species_map = {}
        full_lineage_species_map = {}
        seen_lineages = set()
        if "configuration_table" in self.config_d:
            self.registry.read_configuration_table_species(
                    conf_lineage_species_map=self.config_d["configuration_table"]["lineage_species_map"],
                    conf_constrained_lineages=self.config_d["configuration_table"]["constrained_lineages"],
                    )
        elif SPECIES_LEAFSET_CONSTRAINTS_KEY in self.config_d:
            self.registry.compile_configuration_species_groupings(species_leafset_constraints=self.config_d[SPECIES_LEAFSET_CONSTRAINTS_KEY])

    def write_configuration(
            self,
            output_file,
            output_format,
            output_delimiter="\t",):
        if output_format == "json":
            d = {}
            names = self.registry.preanalysis_constrained_species_lineages_map
            d[SPECIES_LEAFSET_CONSTRAINTS_KEY] = [list(self.registry.preanalysis_constrained_species_lineages_map[n]) for n in names]
            d["species_names"] = list(names)
            if False: #args.output_format == "json-compact":
                json.dump(d, outf)
            else:
                json.dump(d, output_file, indent=4, separators=(',', ': '))
        else:
            output_file.write("{}\n".format(output_delimiter.join(CONFIGURATION_REQUIRED_FIELDS)))
            for lineage_name in self.registry.normalized_tree_lineage_names:
                row = []
                row.append(lineage_name)
                row.append(self.registry.preanalysis_constrained_lineage_species_map.get(lineage_name, "?"))
                if lineage_name in self.registry.preanalysis_constrained_lineage_species_map:
                    row.append("1")
                else:
                    row.append("0")
                output_file.write("{}\n".format(output_delimiter.join(row)))

    def compile_postanalysis_lineage_species_name_map(self, postanalysis_species_leafset_labels):
        return summarize.compile_postanalysis_lineage_species_name_map(
                preanalysis_constrained_lineage_species_map=self.registry.preanalysis_constrained_lineage_species_map,
                postanalysis_species_leafset_labels=postanalysis_species_leafset_labels,
                )

#! /usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import csv
import json
import sys
from delineate import model
from delineate import utility
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

def get_controller(args, name, logger=None):
    controller = Controller(
            name="delineate-estimate",
            is_case_sensitive=getattr(args, "is_case_sensitive", False),
            is_fail_on_extra_tree_lineages=getattr(args, "is_fail_on_extra_tree_lineages", True),
            is_fail_on_extra_configuration_lineages=getattr(args, "is_fail_on_extra_configuration_lineages", True),
            logger=logger)
    controller.read_tree(
            tree_filepath=args.tree_file,
            schema=args.tree_format,
            underflow_protection=getattr(args, "underflow_protection", False))
    controller.parse_configuration_file(
            config_filepath=args.config_file,
            delimiter=None)
    return controller

class Registry(object):

    def __init__(self, **kwargs):
        self.tree_lineages = None
        self.config_file_lineages = None
        self.is_case_sensitive = kwargs.pop("is_case_sensitive", False)
        self.is_fail_on_extra_tree_lineages = kwargs.pop("is_fail_on_extra_tree_lineages", True)
        self.is_fail_on_extra_configuration_lineages = kwargs.pop("is_fail_on_extra_configuration_lineages", True)
        if self.is_case_sensitive:
            self.normalized_lineage_names = {}
            self.normalized_config_lineage_names = {}
        else:
            self.normalized_lineage_names = OrderedCaselessDict()
            self.normalized_config_lineage_names = OrderedCaselessDict()
        self.extra_tree_lineages = []
        self.extra_configuration_lineages = []

    def normalize_lineage_names(self):
        tree_lineage_set = set(self.tree_lineages)
        # tree lineages give the canonical orthography
        for lineage in self.tree_lineages:
            self.normalized_lineage_names[lineage] = lineage
        for lineage in self.config_file_lineages:
            self.normalized_config_lineage_names[lineage] = lineage

    def validate_lineage_names(self, logger):
        for lineage in self.config_file_lineages:
            if lineage not in self.normalized_lineage_names:
                self.extra_configuration_lineages.append(lineage)
        for lineage in self.tree_lineages:
            if lineage not in self.normalized_config_lineage_names:
                self.extra_tree_lineages.append(lineage)

        if self.extra_tree_lineages:
            s1_error_msg = ["{}: {} lineages found on tree but not in configuration data:".format(
                "ERROR" if self.is_fail_on_extra_tree_lineages else "NOTE",
                len(self.extra_tree_lineages))]
            for lidx, label in enumerate(sorted(self.extra_tree_lineages, key=lambda x: x.lower())):
                s1_error_msg.append("    [{: 3d}/{:<3d}] {}".format(lidx+1, len(self.extra_tree_lineages), label))
            s1_error_msg = "\n".join(s1_error_msg)
        else:
            s1_error_msg = ""

        if self.extra_configuration_lineages:
            s2_error_msg = ["{}: {} lineages found in configuration data but not on tree:".format(
                "ERROR" if self.is_fail_on_extra_configuration_lineages else "NOTE",
                len(self.extra_configuration_lineages))]
            for lidx, label in enumerate(sorted(self.extra_configuration_lineages, key=lambda x: x.lower())):
                s2_error_msg.append("    [{: 3d}/{:<3d}] {}".format(lidx+1, len(self.extra_configuration_lineages), label))
            s2_error_msg = "\n".join(s2_error_msg)
        else:
            s2_error_msg = ""

        is_fail = False
        if self.extra_tree_lineages and self.is_fail_on_extra_tree_lineages:
            logger.error(s1_error_msg)
            is_fail = True
        if self.extra_configuration_lineages and self.is_fail_on_extra_configuration_lineages:
            logger.error(s2_error_msg)
            is_fail = True
        if is_fail:
            logger.error("Terminating due to lineage identity conflict error")
            sys.exit(1)

    def compose_report(self):
        msg = []
        msg.append("{} terminal lineages on population tree".format(len(self.tree_lineages)))
        msg.append("{} lineages described in configuration file".format(len(self.config_file_lineages)))

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
        msg = []
        msg.append(("{} fields found in configuration source:".format(len(src_data.fieldnames))))
        for idx, fn in enumerate(src_data.fieldnames):
            msg.append("    [{}/{}] '{}'".format(idx+1, len(src_data.fieldnames), fn))
        self.logger.info("\n".join(msg))
        for required_field in CONFIGURATION_REQUIRED_FIELDS:
            if required_field not in fieldname_set:
                self.logger.error("ERROR: Required field '{}' not found in configuration source".format(required_field))
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
        # self.logger.info("{} lineages in total".format(len(known) + len(unknown)))
        # msg = []
        # msg.append("{} species defined in total".format(len(species_lineage_map)))
        # for spidx, sp in enumerate(species_lineage_map):
        #     msg.append("    [{}/{}] '{}'".format(spidx+1, len(species_lineage_map), sp))
        # msg = "\n".join(msg)
        # self.logger.info(msg)
        # self.logger.info("{} lineages assigned to {} species".format(len(known), len(species_lineage_map)))
        # self.logger.info("{} lineages of unknown species affinities".format(len(unknown)))
        species_leafset_constraints = []
        for key in species_lineage_map:
            species_leafset_constraints.append(species_lineage_map[key])
        assert len(species_leafset_constraints) == len(species_lineage_map)
        self.config_d = {}
        self.config_d[SPECIES_LEAFSET_CONSTRAINTS_KEY] = species_leafset_constraints
        self.config_d["configuration_file"] = {}
        self.config_d["configuration_file"]["lineages"] = known + unknown
        self.config_d["configuration_file"]["constrained_lineages"] = known
        self.config_d["configuration_file"]["unconstrained_lineages"] = unknown
        self.config_d["configuration_file"]["lineage_species_map"] = lineage_species_map
        self.config_d["configuration_file"]["constrained_lineage_species_map"] = {lineage_name:lineage_species_map[lineage_name] for lineage_name in known}
        self.config_d["configuration_file"]["species_lineage_map"] = species_lineage_map
        return self.config_d

    def register_names(self):
        self.register_lineage_names()

    def register_lineage_names(self):
        if self.tree is None:
            raise ValueError("'tree' not set")
        self.registry.tree_lineages = [t.label for t in self.tree.taxon_namespace]
        if "configuration_file" in self.config_d and "lineages" in self.config_d["configuration_file"]:
            self.registry.config_file_lineages = list(self.config_d["configuration_file"]["lineages"])
        else:
            self.registry.config_file_lineages = list(self.tree_lineages)
            # self.registry.config_file_lineages = []
        if not self.registry.config_file_lineages:
            return
        self.registry.normalize_lineage_names()
        self.registry.validate_lineage_names(self.logger)

    def xvalidate_configuration(
            self,
            output_file=None,
            output_format="json",
            is_pretty_print=True,
            output_delimiter="\t",
            is_case_sensitive=False,
            is_fail_on_extra_tree_lineages=False,
            is_fail_on_extra_configuration_lineages=True,
            ):
        msg = []
        self.tree_lineages = None
        config_lineages = None
        if self.tree is not None:
            self.tree_lineages = [t.label for t in self.tree.taxon_namespace]
            msg.append("{} terminal lineages on population self.tree".format(len(self.tree_lineages)))
        else:
            self.tree_lineages = None
        if "configuration_file" in self.config_d and "lineages" in self.config_d["configuration_file"]:
            config_lineages = list(self.config_d["configuration_file"]["lineages"])
            msg.append("{} lineages described in configuration file".format(len(config_lineages)))
        else:
            if self.tree_lineages is not None:
                config_lineages = list(self.tree_lineages)
            else:
                raise ValueError("Tree file required for checking JSON file configuration")
        lineage_case_normalization_map = {}
        if self.tree_lineages is not None and config_lineages is not None:
            self.tree_lineage_set = set(self.tree_lineages)
            config_lineage_set = set(config_lineages)
            if is_case_sensitive:
                self.tree_lineages_lower = {}
                for lineage in self.tree_lineage_set:
                    self.tree_lineages_lower[lineage.lower()] = lineage
                    # lineage_case_normalization_map[lineage] = lineage
                for lineage in config_lineage_set:
                    if lineage.lower() in self.tree_lineages_lower:
                        lineage_case_normalization_map[lineage] = self.tree_lineages_lower[lineage.lower()]
                    else:
                        lineage_case_normalization_map[lineage] = lineage
                config_lineage_set_normalized = set([lineage_case_normalization_map[lineage] for lineage in config_lineage_set])
                config_lineage_set = config_lineage_set_normalized
            s1 = self.tree_lineage_set - config_lineage_set
            if s1:
                # msg.append("NOTE: {} lineages in self.tree not described in configuration file: {}".format(len(s1), ", ".join(s1)))
                s1_error_msg = ["{}: {} lineages found on self.tree but not in configuration file:".format(
                    "ERROR" if is_fail_on_extra_tree_lineages else "NOTE",
                    len(s1))]
                for lidx, label in enumerate(sorted(s1, key=lambda x: x.lower())):
                    s1_error_msg.append("    [{: 3d}/{:<3d}] {}".format(lidx+1, len(s1), label))
                s1_error_msg = "\n".join(s1_error_msg)
            else:
                s1_error_msg = ""
            s2 = config_lineage_set - self.tree_lineage_set
            if s2:
                # msg.append("WARNING: {} lineages in configuration file not found on self.tree: {}".format(len(s2), ", ".join(s2)))
                s2_error_msg = ["{}: {} lineages described in configuration file but not found on tree:".format(
                    "ERROR" if is_fail_on_extra_tree_lineages else "WARNING",
                    len(s2))]
                for lidx, label in enumerate(sorted(s2, key=lambda x: x.lower())):
                    s2_error_msg.append("    [{: 3d}/{:<3d}] {}".format(lidx+1, len(s2), label))
                s2_error_msg = "\n".join(s2_error_msg)
            else:
                s2_error_msg = ""
            if s1 or s2:
                if s1 and self.is_fail_on_extra_tree_lineages:
                    self.logger.error(s1_error_msg)
                elif s1:
                    msg.append(s1_error_msg)
                if s2 and self.is_fail_on_extra_configuration_lineages:
                    self.logger.error(s2_error_msg)
                elif s2:
                    msg.append(s2_error_msg)
                if s1 and self.is_fail_on_extra_tree_lineages:
                    self.logger.error("Exiting due to lineage identity mismatch error (1)")
                    sys.exit(1)
                if s2 and self.is_fail_on_extra_configuration_lineages:
                    self.logger.error("Exiting due to lineage identity mismatch error (2)")
                    sys.exit(1)
            all_lineages = sorted(set(self.tree_lineages + config_lineages))
        elif self.tree_lineages is not None:
            all_lineages = sorted(self.tree_lineages)
        elif config_lineages is not None:
            all_lineages = sorted(config_lineages)

        ### MAPPING

        species_lineage_map = {}
        constrained_lineage_species_map = {}
        full_lineage_species_map = {}
        seen_lineages = set()
        if "configuration_file" in self.config_d:
            # additional info from configurator front-end
            for spp in self.config_d["configuration_file"]["species_lineage_map"]:
                species_lineage_map[spp] = []
                for lineage_name in self.config_d["configuration_file"]["species_lineage_map"][spp]:
                    normalized_lineage_name = lineage_case_normalization_map.get(lineage_name, lineage_name)
                    if normalized_lineage_name in seen_lineages:
                        self.logger.error("ERROR: Duplicate lineage species assignment: '{}'".format(normalized_lineage_name))
                        sys.exit(1)
                    seen_lineages.add(normalized_lineage_name)
                    species_lineage_map[spp].append(normalized_lineage_name)
            constrained_lineage_species_map = {}
            for lineage_name in self.config_d["configuration_file"]["lineage_species_map"]:
                normalized_lineage_name = lineage_case_normalization_map.get(lineage_name, lineage_name)
                full_lineage_species_map[normalized_lineage_name] = self.config_d["configuration_file"]["lineage_species_map"][lineage_name]
                if lineage_name in self.config_d["configuration_file"]["constrained_lineages"]:
                    constrained_lineage_species_map[normalized_lineage_name] = full_lineage_species_map[normalized_lineage_name]
        elif SPECIES_LEAFSET_CONSTRAINTS_KEY in self.config_d:
            species_leafsets = self.config_d[SPECIES_LEAFSET_CONSTRAINTS_KEY]
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
        if self.logger:
            pmsg = "\n".join(["  - {}".format(m) for m in msg])
            self.logger.info("Analysis configuration:\n{}".format(pmsg))
        if output_file is None:
            pass
        elif output_format == "json":
            kwargs = {}
            if is_pretty_print:
                kwargs["sort_keys"] = True
                kwargs["indent"] = 4
            json.dump(self.config_d, output_file, **kwargs)
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
                # "messages": msg,
            }

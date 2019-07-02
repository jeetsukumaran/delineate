#! /usr/bin/env python
# -*- coding: utf-8 -*-

def compile_postanalysis_lineage_species_name_map(
        preanalysis_constrained_lineage_species_map,
        postanalysis_species_leafset_labels,
        is_validate_species_group_consistency=True,
        ):
    lnsp_map = {}
    new_sp_idx = 0
    existing_sp_names = set(preanalysis_constrained_lineage_species_map.values())
    processed_sp_names = set()
    for leafset in postanalysis_species_leafset_labels:
        sp_id = None
        for lineage in leafset:
            try:
                assigned_sp_id = preanalysis_constrained_lineage_species_map[lineage]
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


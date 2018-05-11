#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from collections import defaultdict
import copy
import random
import math
import dendropy

################################################################################
## Functions in support of calculating the joint probability

def _gen_from_sub_list(sub_list, open_el):
    # print('sub_list={}'.format(sub_list))
    sub_list.sort()
    sub_list = tuple(sub_list)
    new_ot_ind = -1 if open_el is None else sub_list.index(open_el)
    return _Partition(data=(new_ot_ind, sub_list))

def _del_part_maps(nd):
    try:
        delattr(nd, 'tipward_part_map')
    except AttributeError:
        pass
    try:
        delattr(nd, 'rootward_part_map')
    except AttributeError:
        pass

def _merge_into_first(src_and_dest_dict, second):
    first = copy.copy(src_and_dest_dict)
    src_and_dest_dict.clear()
    dest = src_and_dest_dict
    for fpart, fprob in first.items():
        # print('fpart = {}'.format(fpart))
        for spart, sprob in second.items():
            # print('spart = {}'.format(spart))
            dest[fpart.create_extension(spart)] += fprob * sprob
    return dest

def _create_closed_map(map_with_open, remain_open_prob, allow_closing=True):
    ret = defaultdict(float)
    if allow_closing:
        close_prob = 1.0 - remain_open_prob
        for part, prob in map_with_open.items():
            if part.is_open:
                ret[part] += prob * remain_open_prob
                ret[part.create_closed()] += prob * close_prob
            else:
                ret[part] += prob
    else:
        for part, prob in map_with_open.items():
            if part.is_open:
                ret[part] += prob * remain_open_prob
    return ret

def _enforce_constraints(partition_table, constraints):
    consp_constraints = constraints.get('conspecific', [])
    not_consp_constraints = constraints.get('not_conspecific', [])
    to_del = []
    for k in partition_table.keys():
        if k.violates_constraints(consp_constraints, not_consp_constraints):
            to_del.append(k)
    for k in to_del:
        del partition_table[k]
    return len(to_del)


# noinspection PyProtectedMember
class _Partition(object):
    __slots__ = '_data'

    @staticmethod
    def compile_lookup_key(group_of_groups):
        return frozenset( frozenset(i) for i in group_of_groups )

    def __init__(self, leaf_label=None, data=None):
        if data is not None:
            self._data = data
        else:
            assert leaf_label is not None
            _label_subsets = [[leaf_label]]
            open_index = 0
            as_l = [list(i) for i in _label_subsets]
            for i in as_l:
                i.sort()
            as_l.sort()
            d = [open_index, tuple([tuple(i) for i in as_l])]
            self._data = tuple(d)

    @property
    def index_of_open_el(self):
        return self._data[0]

    @property
    def _label_subsets(self):
        return self._data[1]

    def violates_constraints(self, consp_constraints, not_consp_constraints):
        s_t = self._label_subsets
        for cc in consp_constraints:
            first, second = cc
            ff, sf = False, False
            for sub in s_t:
                if first in sub:
                    ff = True
                    if second in sub:
                        break
                    elif sf:
                        return True
                    if sf:
                        break
                if second in sub:
                    sf = True
                    if ff:
                        return True
        for cc in not_consp_constraints:
            first, second = cc
            for sub in s_t:
                if first in sub:
                    if second in sub:
                        return True
                    break
                elif second in sub:
                    break
        return False

    def create_closed(self):
        asl = list(self._data)
        asl[0] = -1
        return _Partition(data=tuple(asl))

    def create_extension(self, other):
        new_sub = []
        new_ot = None
        s_subs = self._label_subsets
        o_subs = other._label_subsets
        if self.is_open:
            s_opind = self.index_of_open_el
            s_op_el = s_subs[s_opind]
            if other.is_open:
                for el in s_subs:
                    if el is not s_op_el:
                        new_sub.append(el)
                o_opind = other.index_of_open_el
                o_op_el = o_subs[o_opind]
                for el in o_subs:
                    if el is not o_op_el:
                        new_sub.append(el)
                ls_op_el = list(s_op_el)
                ls_op_el.extend(o_op_el)
                ls_op_el.sort()
                new_ot = tuple(ls_op_el)
                new_sub.append(new_ot)
                return _gen_from_sub_list(new_sub, new_ot)
            new_ot = s_op_el
        elif other.is_open:
            return other.create_extension(self)
        new_sub.extend(s_subs)
        new_sub.extend(o_subs)
        return _gen_from_sub_list(new_sub, new_ot)

    @property
    def is_open(self):
        return self.index_of_open_el >= 0

    def __eq__(self, other):
        return self._data == other._data

    def __hash__(self):
        return hash(self._data)

    def __str__(self):
        as_str_l = [str(i) for i in self._label_subsets]
        if self.is_open:
            as_str_l[self.index_of_open_el] = as_str_l[self.index_of_open_el] + '*'
        return '_Partition[{}]'.format(' | '.join(as_str_l))

    @property
    def set_notation(self):
        as_l = [','.join(i) for i in self._label_subsets]
        return ''.join(['{', '}{'.join(as_l), '}'])

    def lookup_key(self):
        return _Partition.compile_lookup_key(self._data[1])

################################################################################
## _Cache class

class _Cache(object):

    def __init__(self):
        self._calc_fns = {}
        self._values = {}

    def clear(self, key=None):
        if key is None:
            self._values = {}
        else:
            self._values.pop(key, None)

    def recalc(self, key=None):
        if key is None:
            self._values = {}
            for k in self._calc_fns:
                self._values[k] = self._calc_fns[k]()
        else:
            self._values[key] = self._calc_fns[key]

    def add(self, key, calc_fn):
        self._calc_fns[key] = calc_fn
        self._values.pop(key, None)

    def __getitem__(self, key):
        try:
            return self._values[key]
        except KeyError:
            self._values[key] = self._calc_fns[key]()
            return self._values[key]

################################################################################
## Enum class

class SF(object):
    """
    Enum class to name bits associated with being anc of none, some and all
    selected tips.
    """
    UNSET   = 0   # has no descendants that are selected
    SEL_DES = 1   # anc to some selected tips
    CA_BIT  = 2   # bit representing "is a common anc to all of the selected tips"
    CA_FLAG = 3   # flag for a common ancestor of all selected tips


################################################################################
## LineageEdge

class LineageEdge(dendropy.Edge):

    def __init__(self, tree, **kwargs):
        dendropy.Edge.__init__(self, **kwargs)
        self.tree = tree

    def _get_length(self):
        return self._length
    def _set_length(self, length):
        self._length = length
        # to do! self.tree.invalidate_cache(self)
    length = property(_get_length, _set_length)

################################################################################
## LineageNode

class LineageNode(dendropy.Node):

    def __init__(self, tree, **kwargs):
        self.tree = tree
        self.marginal_prob_calc = {}
        dendropy.Node.__init__(self, **kwargs)

    def edge_factory(self, **kwargs):
        return LineageEdge(tree=self.tree, **kwargs)

################################################################################
## LineageTree

class LineageTree(dendropy.Tree):

    ################################################################################
    ## Lifecycle and Structure

    def __init__(self, *args, **kwargs):
        self._speciation_completion_rate = None
        self._setup_cache()
        dendropy.Tree.__init__(self, *args, **kwargs)

    def node_factory(self, *args, **kwargs):
        return LineageNode(tree=self, **kwargs)

    ################################################################################
    ## Properties

    def _get_speciation_completion_rate(self):
        return self._speciation_completion_rate

    def _set_speciation_completion_rate(self, value):
        self._speciation_completion_rate = value
        self.invalidate_cache(self)

    speciation_completion_rate = property(_get_speciation_completion_rate, _set_speciation_completion_rate)

    ################################################################################
    ## Cache

    def _setup_cache(self):
        self._joint_probability_cache = _Cache()
        self._joint_probability_cache.add(
                "label_partition_probability_map",
                self.calc_label_partition_probability_map)

    def invalidate_cache(self, o):
        # o = object that changed that required cache invalidation
        self._joint_probability_cache.clear()

    ################################################################################
    ## Marginal Probability

    def calc_marginal_probability_of_species(self, selected_tip_labels):
        """
        Calculates the marginal probability that there is a "good" species with the tip labels
        that correspond to the set `selected_tip_labels`.
        """
        num_sel = len(selected_tip_labels)
        sel_as_flag = SF.CA_FLAG if num_sel == 1 else SF.SEL_DES
        total_prob = 0.0
        for nd in self.postorder_node_iter():
            if nd.is_leaf():
                if nd.taxon.label in selected_tip_labels:
                    nd.marginal_prob_calc["num_sel"] = 1
                    nd.marginal_prob_calc["anc_status"] = sel_as_flag
                    nd.marginal_prob_calc["accum_prob"] = 1.0
                else:
                    nd.marginal_prob_calc["num_sel"] = 0
                    nd.marginal_prob_calc["anc_status"] = SF.UNSET
                    nd.marginal_prob_calc["accum_prob"] = 0.0
            else:
                nd.marginal_prob_calc["num_sel"] = 0
                for c in nd.child_nodes():
                    nd.marginal_prob_calc["num_sel"] += c.marginal_prob_calc["num_sel"]
                if nd.marginal_prob_calc["num_sel"] == 0:
                    nd.marginal_prob_calc["anc_status"] = SF.UNSET
                elif nd.marginal_prob_calc["num_sel"] == num_sel:
                    nd.marginal_prob_calc["anc_status"] = SF.CA_FLAG
                else:
                    nd.marginal_prob_calc["anc_status"] = SF.SEL_DES
                total_prob += self._marginal_species_prob_accum_prob(nd, self._speciation_completion_rate)
        total_prob += self.seed_node.marginal_prob_calc["accum_prob"]
        return total_prob

    def _marginal_species_prob_accum_prob(self, nd, good_sp_rate):
        """
        Fills in the accum_prob slot for nd, and returns any contribution to
        the probability of the selected taxa being a good species.
        """
        ap = 1.0
        ret = 0.0
        for child in nd.child_nodes():
            c_brlen = child.edge.length
            scaled_brlen = c_brlen * good_sp_rate
            prob_no_sp = math.exp(-scaled_brlen)
            prob_sp = 1.0 - prob_no_sp
            if child.marginal_prob_calc["anc_status"] & SF.SEL_DES:
                if child.marginal_prob_calc["anc_status"] & SF.CA_BIT:
                    ret = prob_sp * child.marginal_prob_calc["accum_prob"]
                contrib = prob_no_sp * child.marginal_prob_calc["accum_prob"]
            else:
                contrib = prob_sp + prob_no_sp * child.marginal_prob_calc["accum_prob"]
            ap *= contrib
        nd.marginal_prob_calc["accum_prob"] = ap
        return ret

    ################################################################################
    ## Joint Probability

    def calc_joint_probability_of_species(self, taxon_labels):
        if not isinstance(taxon_labels, set):
            taxon_labels = frozenset(frozenset(i) for i in taxon_labels)
        prob = self._joint_probability_cache["label_partition_probability_map"][taxon_labels]
        return prob

    def calc_label_partition_probability_map(self):
        partition_probability_map = self._calc_all_joint_sp_probs(good_sp_rate=self._speciation_completion_rate)
        return partition_probability_map

    def _calc_all_joint_sp_probs(self, good_sp_rate):
        for nd in self.postorder_node_iter():
            if nd.is_leaf():
                nd.tipward_part_map = defaultdict(float)
                nd.tipward_part_map[_Partition(leaf_label=nd.taxon.label)] = 1.0
            else:
                children = nd.child_nodes()
                nd.tipward_part_map = children[0].rootward_part_map
                _del_part_maps(children[0])
                constraints = getattr(nd, 'sp_constraints', None)
                for c in children[1:]:
                    _merge_into_first(nd.tipward_part_map, c.rootward_part_map)
                    _del_part_maps(c)
                    if constraints is not None:
                        _enforce_constraints(nd.tipward_part_map, constraints)
            if nd is self.seed_node:
                break
            c_brlen = nd.edge.length
            scaled_brlen = c_brlen * good_sp_rate
            prob_no_sp = math.exp(-scaled_brlen)
            nd.rootward_part_map = _create_closed_map(nd.tipward_part_map,
                                                      prob_no_sp,
                                                      getattr(nd, 'speciation_allowed', True))
        final_part_map = defaultdict(float)
        if False:
            # use _Partition as key
            for part, prob in self.seed_node.tipward_part_map.items():
                k = part.create_closed() if part.is_open else part
                final_part_map[k] += prob
        else:
            # use lookup key as key
            for part, prob in self.seed_node.tipward_part_map.items():
                k = part.create_closed() if part.is_open else part
                final_part_map[k.lookup_key()] += prob
        _del_part_maps(self.seed_node)
        return final_part_map


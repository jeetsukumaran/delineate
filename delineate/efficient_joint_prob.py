#!/usr/bin/env python
from __future__ import print_function
from collections import defaultdict
import copy
import dendropy
import math


# noinspection PyProtectedMember
class Partition(object):
    __slots__ = '_data'

    def __init__(self, leaf_label=None, data=None):
        if data is not None:
            self._data = data
        else:
            assert leaf_label is not None
            subsets = [[leaf_label]]
            open_index = 0
            as_l = [list(i) for i in subsets]
            for i in as_l:
                i.sort()
            as_l.sort()
            d = [open_index, tuple([tuple(i) for i in as_l])]
            self._data = tuple(d)

    @property
    def index_of_open_el(self):
        return self._data[0]

    @property
    def subsets(self):
        return self._data[1]

    def create_closed(self):
        asl = list(self._data)
        asl[0] = -1
        return Partition(data=tuple(asl))

    def create_extension(self, other):
        new_sub = []
        new_ot = None
        s_subs = self.subsets
        o_subs = other.subsets
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
        as_str_l = [str(i) for i in self.subsets]
        if self.is_open:
            as_str_l[self.index_of_open_el] = as_str_l[self.index_of_open_el] + '*'
        return 'Partition[{}]'.format(' | '.join(as_str_l))

    @property
    def set_notation(self):
        as_l = [','.join(i) for i in self.subsets]
        return ''.join(['{', '}{'.join(as_l), '}'])


def _gen_from_sub_list(sub_list, open_el):
    # print('sub_list={}'.format(sub_list))
    sub_list.sort()
    sub_list = tuple(sub_list)
    new_ot_ind = -1 if open_el is None else sub_list.index(open_el)
    return Partition(data=(new_ot_ind, sub_list))


def del_part_maps(nd):
    delattr(nd, 'tipward_part_map')
    delattr(nd, 'rootward_part_map')


def merge_into_first(src_and_dest_dict, second):
    first = copy.copy(src_and_dest_dict)
    src_and_dest_dict.clear()
    dest = src_and_dest_dict
    for fpart, fprob in first.items():
        # print('fpart = {}'.format(fpart))
        for spart, sprob in second.items():
            # print('spart = {}'.format(spart))
            dest[fpart.create_extension(spart)] += fprob * sprob
    return dest


def create_closed_map(map_with_open, remain_open_prob):
    ret = defaultdict(float)
    close_prob = 1.0 - remain_open_prob
    for part, prob in map_with_open.items():
        if part.is_open:
            ret[part] += prob * remain_open_prob
            ret[part.create_closed()] += prob * close_prob
        else:
            ret[part] += prob
    return ret


def calc_all_joint_sp_probs(tree, good_sp_rate):
    """Calculates the marginal probability that there is a "good" species with the tip labels
    that correspond to the set `selected_tip_labels`.
    """
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            nd.tipward_part_map = defaultdict(float)
            nd.tipward_part_map[Partition(leaf_label=nd.taxon.label)] = 1.0
        else:
            children = nd.child_nodes()
            nd.tipward_part_map = children[0].rootward_part_map
            del_part_maps(children[0])
            for c in children[1:]:
                merge_into_first(nd.tipward_part_map, c.rootward_part_map)
                del_part_maps(c)
        if nd is tree.seed_node:
            break
        c_brlen = nd.edge.length
        scaled_brlen = c_brlen * good_sp_rate
        prob_no_sp = math.exp(-scaled_brlen)
        nd.rootward_part_map = create_closed_map(nd.tipward_part_map, prob_no_sp)
    final_part_map = defaultdict(float)
    for part, prob in tree.seed_node.tipward_part_map.items():
        k = part.create_closed() if part.is_open else part
        final_part_map[k] += prob
    return final_part_map


def main(tree_filename, good_sp_rate):
    tree = dendropy.Tree.get(path=tree_filename, schema="newick")
    r = calc_all_joint_sp_probs(tree, good_sp_rate)
    for part, prob in r.items():
        print('Pr( {} ) = {}'.format(part.set_notation, prob))


if __name__ == '__main__':
    import sys

    try:
        rate = 1.0
        filename = sys.argv[1]
        if len(sys.argv) > 2:
            assert len(sys.argv) == 3
            rate = float(sys.argv[2])
        assert rate > 0.0
    except:
        sys.exit('''Expecting up to 2 args: 
    1. the filepath to a rooted newick tree with branch lengths, and
    2. (optionally) a rate of good speciation events (branch length multiplier), and
''')
    main(filename, rate)

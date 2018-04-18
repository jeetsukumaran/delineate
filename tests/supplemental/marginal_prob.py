#!/usr/bin/env python
from __future__ import print_function
import dendropy
import math


class SF(object):
    """Enum class to name bits associated with being anc of none, some and all selected tips."""
    UNSET = 0    # has no descendants that are selected
    SEL_DES = 1  # anc to some selected tips
    CA_BIT = 2   # bit representing "is a common anc to all of the selected tips"
    CA_FLAG = 3  # flag for a common ancestor of all selected tips


def calc_marginal_probability_of_species(tree, selected_tip_labels, good_sp_rate):
    """Calculates the marginal probability that there is a "good" species with the tip labels
    that correspond to the set `selected_tip_labels`.
    """
    num_sel = len(selected_tip_labels)
    sel_as_flag = SF.CA_FLAG if num_sel == 1 else SF.SEL_DES
    total_prob = 0.0
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            if nd.taxon.label in selected_tip_labels:
                nd.num_sel = 1
                nd.anc_status = sel_as_flag
                nd.accum_prob = 1.0
            else:
                nd.num_sel = 0
                nd.anc_status = SF.UNSET
                nd.accum_prob = 0.0
        else:
            nd.num_sel = 0
            for c in nd.child_nodes():
                nd.num_sel += c.num_sel
            if nd.num_sel == 0:
                nd.anc_status = SF.UNSET
            elif nd.num_sel == num_sel:
                nd.anc_status = SF.CA_FLAG
            else:
                nd.anc_status = SF.SEL_DES
            total_prob += accum_prob(nd, good_sp_rate)
    total_prob += tree.seed_node.accum_prob
    return total_prob


def accum_prob(nd, good_sp_rate):
    """Fills in the accum_prob slot for nd, and returns any contribution to the probability of
    the selected taxa being a good species.
    """
    ap = 1.0
    ret = 0.0
    for child in nd.child_nodes():
        c_brlen = child.edge.length
        scaled_brlen = c_brlen * good_sp_rate
        prob_no_sp = math.exp(-scaled_brlen)
        prob_sp = 1.0 - prob_no_sp
        if child.anc_status & SF.SEL_DES:
            if child.anc_status & SF.CA_BIT:
                ret = prob_sp * child.accum_prob
            contrib = prob_no_sp * child.accum_prob
        else:
            contrib = prob_sp + prob_no_sp * child.accum_prob
        ap *= contrib
    nd.accum_prob = ap
    return ret


def main(tree_filename, good_sp_rate, selected_tip_labels):
    tree = dendropy.Tree.get(path=tree_filename, schema="newick")
    selected_nodes = []
    labels_found = set()
    for tip in tree.leaf_node_iter():
        tl = tip.taxon.label
        if tl in selected_tip_labels:
            if tl in labels_found:
                sys.exit('Tip label "{}" occurred twice in the tree!\n'.format(tl))
            labels_found.add(tl)
            selected_nodes.append(tip)
    if labels_found != selected_tip_labels:
        d = selected_tip_labels - labels_found
        sys.exit('Not all tips were found. Missing: "{}"\n'.format('", "'.join(list(d))))
    prob_good = calc_marginal_probability_of_species(tree, selected_tip_labels, good_sp_rate)
    stn = list(selected_tip_labels)
    stn.sort()
    stl = ','.join(stn)
    print('Pr({' + stl + '}) = ' + str(prob_good))


if __name__ == '__main__':
    import sys

    try:
        rate = 1.0
        filename = sys.argv[1]
        assert len(sys.argv) > 3
        rate = float(sys.argv[2])
        assert rate > 0.0
        selected_tip_label_list = sys.argv[3:]
        selected_tip_label_set = set(selected_tip_label_list)
        if len(selected_tip_label_set) < len(selected_tip_label_list):
            sys.stderr.write('WARN: some tip labels were repeated in the command line.\n')
    except:
        sys.exit('''Expecting up to at least 3 args:
    1. the filepath to a rooted newick tree with branch lengths,
    2. a rate of good speciation events (branch length multiplier), and
    3. a series of taxa labels that designate the possible conspecific lineages.
''')
    main(filename, rate, selected_tip_label_set)

#!/usr/bin/env python
from __future__ import print_function
import dendropy
import math


def format_partition(p):
    return ''.join(['{', '}{'.join([','.join(sub) for sub in p]), '}'])


# Note: node flagging somewhat different from notes to unite mono and non-mono:
#  anc_status flag is a union of bits which mean:
#   1 = Some Des are Selected nodes
#   2 = Some Des are Unselected nodes
#   4 = All Selected nodes are descendants.
# So in the monophyletic case the statuses will be: 1, 2, 5 (selected MRCA) and 7 (anc of MRCA)
# While in the non-monophyletic form, so statuses of 3 will also be encountered
# So the meaning of entire flags:
#   1 (SEL_DES) => anc of some (but not all) selected nodes
#   2 (UNSEL_DES) => anc of only unselected nodes
#   3 (MIXED_DES) => anc of some (but not all) selected nodes and some unselected nodes
#   4 (MONO_MRCA_BIT) IMPOSSIBLE as flag value
#   5 (SEL_MRCA) => MRCA of all of the selected nodes and no unselected nodes.
#   6 IMPOSSIBLE as a flag value
#   7 (COMM_ANC_MIXED) => ancestor (not necessarily MR anc) of all selected nodes and some unselected
#       if the selected nodes are monophyletic, this is an ancestor of SEL_MRCA. If they are not
#       monophyletic, this could be the MRCA of the selected node or it ancestor.
class SF(object):
    # Flags. 6 is impossible.
    SEL_DES = 1  # has some selected des
    UNSEL_DES = 2
    MIXED_DES = 3
    MONO_MRCA_BIT = 4
    SEL_MRCA = 5
    COMM_ANC_MIXED = 7


def calc_prob_good_species(tree, selected_tips, good_sp_rate):
    # Add an attribute to each tip for the number of tips "selected" and "unselected" at
    #   or below this node.
    for tip in tree.leaf_node_iter():
        tip.num_sel = 0
        tip.num_unsel = 1
        tip.anc_status = SF.UNSEL_DES
        tip.accum_prob = 0.0
    one_tip_sel = len(selected_tips) == 1
    sel_as_flag = SF.SEL_MRCA if one_tip_sel else SF.SEL_DES
    for tip in selected_tips:
        tip.num_sel = 1
        tip.num_unsel = 0
        tip.anc_status = sel_as_flag
        tip.accum_prob = 1.0
    num_sel = len(selected_tips)
    is_mono = False
    total_prob = 0.0
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            continue
        nd.num_sel, nd.num_unsel = 0, 0
        for c in nd.child_nodes():
            nd.num_unsel += c.num_unsel
            nd.num_sel += c.num_sel
        if nd.num_sel == 0:
            assert nd.num_unsel > 0
            nd.anc_status = SF.UNSEL_DES
        elif nd.num_unsel == 0:
            if nd.num_sel == num_sel:
                is_mono = True
                nd.anc_status = SF.SEL_MRCA
            else:
                nd.anc_status = SF.SEL_DES
        elif nd.num_sel == num_sel:
            nd.anc_status = SF.COMM_ANC_MIXED
        else:
            assert not is_mono
            nd.anc_status = SF.MIXED_DES
        total_prob += accum_prob(nd, good_sp_rate)
    total_prob += tree.seed_node.accum_prob
    return total_prob


def accum_prob(nd, good_sp_rate):
    # Calculate prob...
    ap = 1.0
    ret = 0.0
    cl = []
    for child in nd.child_nodes():
        if not child.label:
            child.label = child.taxon.label
        cl.append(child.label)
        c_brlen = child.edge.length
        scaled_brlen = c_brlen * good_sp_rate
        prob_no_sp = math.exp(-scaled_brlen)
        prob_sp = 1.0 - prob_no_sp
        if child.anc_status & SF.SEL_DES:
            if child.anc_status & SF.MONO_MRCA_BIT:
                ret = prob_sp * child.accum_prob
            contrib = prob_no_sp * child.accum_prob
        else:
            contrib = prob_sp + prob_no_sp * child.accum_prob
        ap *= contrib
        # fmt = 'child {} flag= {} brlen = {} prob_sp = {} child.accum_prob ={} par.accum_prob = {} ret={}'
        # print(fmt.format(child.label, child.anc_status, c_brlen, prob_sp, child.accum_prob, ap, ret))

    nd.label = '({})'.format(','.join(cl))
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
    prob_good = calc_prob_good_species(tree, selected_nodes, good_sp_rate)
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

#!/usr/bin/env python
from __future__ import print_function
import dendropy
import math

def format_partition(p):
    return ''.join(['{', '}{'.join([','.join(sub) for sub in p]), '}'])

def exhaustive_good_speciation_events(tree, rate):
    node_list = []
    prob_no_sp = []
    prob_sp = []
    leaf_list = []
    for nd in tree.preorder_node_iter():
        if nd is tree.seed_node:
            continue
        scaled_rate = nd.edge.length * rate
        assert scaled_rate >= 0.0
        pns = math.exp(-scaled_rate)
        prob_no_sp.append(pns)
        prob_sp.append(1.0 - pns)
        node_list.append(nd)
        if nd.is_leaf():
            leaf_list.append(nd)
    max_bits = 1 << len(prob_no_sp)
    scenario_index = 0
    root = tree.seed_node
    joint_consp_dict = {}
    while scenario_index < max_bits:
        root.sp_set = set()
        scenario_prob = 1.0
        curr_bit = 1
        for nd, pns, ps in zip(node_list, prob_no_sp, prob_sp):
            if curr_bit & scenario_index:
                p = ps
                nd.sp_set = set()
            else:
                p = pns
                nd.sp_set = nd.parent_node.sp_set
            if nd.is_leaf():

                nd.sp_set.add(nd.taxon.label)
            scenario_prob *= p
            #print(nd, pns, ps, curr_bit, scenario_index, scenario_prob)
            curr_bit <<= 1
        consp_list = []
        consp_set = set()
        for nd in leaf_list:
            sl = list(nd.sp_set)
            sl.sort()
            st = tuple(sl)
            if st not in consp_set:
                consp_set.add(st)
                consp_list.append(st)
        consp_list.sort()
        ct = tuple(consp_list)
        prev = joint_consp_dict.get(ct, 0.0)
        joint_consp_dict[ct] = prev + scenario_prob
        scenario_index += 1
    jcdk = [(len(k), k) for k in joint_consp_dict.keys()]
    jcdk.sort()
    print('Joint statements:')
    for ljc, jc in jcdk:
        print('Pr( {} ) = {}'.format(format_partition(jc), joint_consp_dict[jc]))
    marg_to_prob = {}
    for jc, pr in joint_consp_dict.items():
        for el in jc:
            prev = marg_to_prob.get(el, 0)
            marg_to_prob[el] = prev + pr
    mk = [(len(k), k) for k in marg_to_prob.keys()]
    mk.sort()
    print('marginal statements:')
    for lk, k in mk:
        print('Pr({' + ','.join(k) + '}) = ' + str(marg_to_prob[k]))

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
        sys.exit('Expecting up to 2 args: a newick file name and a rate of good speciation events.')
    tree = dendropy.Tree.get(path=filename, schema="newick")
    exhaustive_good_speciation_events(tree, rate)

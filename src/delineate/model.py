#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import random
import math
import dendropy

class SF(object):
    """
    Enum class to name bits associated with being anc of none, some and all
    selected tips.
    """
    UNSET   = 0   # has no descendants that are selected
    SEL_DES = 1   # anc to some selected tips
    CA_BIT  = 2   # bit representing "is a common anc to all of the selected tips"
    CA_FLAG = 3   # flag for a common ancestor of all selected tips

class LineageTree(dendropy.Tree):

    def __init__(self, *args, **kwargs):
        dendropy.Tree.__init__(self, *args, **kwargs)

    def calc_prob_good_species(self, selected_tip_labels, good_sp_rate):
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
                total_prob += self._marginal_species_prob_accum_prob(nd, good_sp_rate)
        total_prob += self.seed_node.accum_prob
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
            if child.anc_status & SF.SEL_DES:
                if child.anc_status & SF.CA_BIT:
                    ret = prob_sp * child.accum_prob
                contrib = prob_no_sp * child.accum_prob
            else:
                contrib = prob_sp + prob_no_sp * child.accum_prob
            ap *= contrib
        nd.accum_prob = ap
        return ret

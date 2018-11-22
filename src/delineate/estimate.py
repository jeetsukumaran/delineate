#! /usr/bin/env python
# -*- coding: utf-8 -*-

import math
import sys
try:
    import scipy.optimize
except ImportError:
    pass

class SpeciationCompletionRateMaximumLikelihoodEstimator(object):

    def __init__(self,
            tree,
            species_leafset_labels,
            initial_speciation_rate,
            min_speciation_rate,
            max_speciation_rate):
        self.tree = tree
        self.species_leafset_labels = species_leafset_labels
        self.initial_speciation_rate = initial_speciation_rate
        self.min_speciation_rate = min_speciation_rate
        self.max_speciation_rate = max_speciation_rate
        assert self.min_speciation_rate > 0.0
        assert self.min_speciation_rate <= self.max_speciation_rate
        assert self.min_speciation_rate <= self.initial_speciation_rate
        assert self.max_speciation_rate >= self.initial_speciation_rate
        self.tree.set_node_constraints(species_leafset_labels=self.species_leafset_labels)

    def _estimate(self,
            f,
            initial_val,
            min_val,
            max_val):
        assert min_val > 0.0
        assert min_val <= max_val
        assert min_val <= initial_val
        assert max_val >= initial_val
        est_result = scipy.optimize.minimize_scalar(f, method="bounded", bounds=(min_val, max_val))
        return est_result.x, est_result.fun
        # brac_res = scipy.optimize.bracket(f,
        #         xa=min_val,
        #         xb=max_val,
        #         )
        # b = brac_res[:3]
        # if b[0] <= 0:
        #     b[0] = 1e-8
        # sys.stderr.write("brackets: {}\n".format(b))
        # while True:
        #     try:
        #         # est_result = scipy.optimize.brent(f, brack=b, full_output=True)
        #         est_result = scipy.optimize.minimize_scalar(f, bounds=b)
        #         break
        #     except ValueError:
        #         # weird bracket interval; default to min/max bounds
        #         b = (min_val, initial_val, max_val)
        #         est_result = scipy.optimize.brent(f, brack=b, full_output=True)
        # value_estimate = est_result[0]
        # value_estimate_prob = est_result[1]
        # return value_estimate, value_estimate_prob

    def estimate_speciation_rate(self):
        if len(self.species_leafset_labels) == 1:
            speciation_completion_rate_estimate = 0.0
            self.tree.speciation_completion_rate = speciation_completion_rate_estimate
            speciation_completion_rate_estimate_prob = self.tree.calc_joint_probability_of_species(species_leafset_labels=self.species_leafset_labels)
        elif self.tree.all_monotypic:
            speciation_completion_rate_estimate = float('inf')
            self.tree.speciation_completion_rate = speciation_completion_rate_estimate
            speciation_completion_rate_estimate_prob = self.tree.calc_joint_probability_of_species(species_leafset_labels=self.species_leafset_labels)
        else:
            def f(x, *args):
                self.tree.speciation_completion_rate = x
                return -1 * self.tree.calc_joint_probability_of_species(species_leafset_labels=self.species_leafset_labels)
            x1, x2 = self._estimate(f=f,
                    initial_val=self.initial_speciation_rate,
                    min_val=self.min_speciation_rate,
                    max_val=self.max_speciation_rate,
                    )
            speciation_completion_rate_estimate = x1
            speciation_completion_rate_estimate_prob = -1 * x2
            return speciation_completion_rate_estimate, math.log(speciation_completion_rate_estimate_prob)

    def estimate_confidence_interval(self, mle_speciation_rate, max_lnl):
        def f0(x, *args):
            self.tree.speciation_completion_rate = x
            prob = self.tree.calc_joint_probability_of_species(species_leafset_labels=self.species_leafset_labels)
            try:
                return abs(max_lnl - 1.96 - math.log(prob))
            except ValueError:
                sys.stderr.write("x = {}, prob = {}\n".format(x, prob))
                raise
        min_val = self.min_speciation_rate
        max_val = mle_speciation_rate - 1e-8
        initial_val = min_val + ((max_val - min_val) / 2.0)
        ci_low, _ = self._estimate(f0,
                initial_val=initial_val,
                min_val=min_val,
                max_val=max_val)
        min_val = mle_speciation_rate + 1e-8
        max_val = self.max_speciation_rate
        initial_val = min_val + ((max_val - min_val) / 2.0)
        ci_high, _ = self._estimate(f0,
                initial_val=initial_val,
                min_val=min_val,
                max_val=max_val)
        return ci_low, ci_high


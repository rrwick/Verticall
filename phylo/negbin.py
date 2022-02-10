"""
Copyright 2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/XXXXXXXXX

This file is part of XXXXXXXXX. XXXXXXXXX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. XXXXXXXXX is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with XXXXXXXXX.
If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
from scipy.stats import nbinom

from .misc import get_mean, get_variance


def get_n_and_p(mean, variance):
    """
    Converts the intuitive mean and variance metrics to the less-intuitive n and p parameters
    needed to define a negative binomial distribution.
    """
    p = mean / variance
    n = mean * p / (1.0 - p)
    return n, p


def fit_negbin_to_distribution(masses):
    """
    This function takes an empirical discrete distribution and tries to find a negative binomial
    distribution which best fits the data.

    In addition to n and p (the negative binomial parameters), this function also adjusts the
    vertical scale (from 0.5 to 1.0) because the empirical distribution may be multimodal, and we
    want to fit the negative binomial distribution to only the strongest mode, which we assume is
    at least half of the distribution.

    Empirical points below the negative binomial are penalised more than empirical points above the
    negative binomial. This is because in a multimodal distribution we expect other peaks to rise
    above the strongest mode.
    """
    masses = np.array(masses)
    neg_penalty = 2.0  # Controls how much being low is worse than being high.

    # We start with n/p based on the mean/variance of the distribution and a vertical scale
    # of one. For unimodal distributions, this should get use pretty close to the correct answer.
    mean, variance = get_mean(masses), get_variance(masses)
    if variance <= mean:
        variance = mean * 1.1
    n, p = get_n_and_p(mean, variance)
    vert = 1.0
    best_n, best_p, best_vert = n, p, vert
    best_fit = assess_fit(masses, n, p, vert, neg_penalty=neg_penalty)

    n_step, p_step, vert_step = 0.1, 0.1, 0.1
    while True:
        improvement = False
        mutated_parameters = get_mutated_parameters(best_n, best_p, best_vert,
                                                    n_step, p_step, vert_step)
        for n, p, vert in mutated_parameters:
            fit = assess_fit(masses, n, p, vert, neg_penalty=neg_penalty)
            if best_fit - fit > 0.0000001:
                best_n, best_p, best_vert = n, p, vert
                best_fit = fit
                improvement = True
        if not improvement:
            n_step /= 2.0
            p_step /= 2.0
            vert_step /= 2.0
        if n_step < 0.0000001:
            break

    return best_n, best_p, best_vert


def get_mutated_parameters(n, p, vert, n_step, p_step, vert_step):
    mutated_parameters = []
    if n_step > 0.0:
        high_n = n + n_step
        low_n = n - n_step
        mutated_parameters.append((high_n, p, vert))
        if low_n > 0.0:
            mutated_parameters.append((low_n, p, vert))
    if p_step > 0.0:
        high_p = p + p_step
        low_p = p - p_step
        mutated_parameters.append((n, high_p, vert))
        if low_p > 0.0:
            mutated_parameters.append((n, low_p, vert))
    if vert_step > 0.0:
        high_vert = vert + vert_step
        low_vert = vert - vert_step
        if high_vert <= 1.0:
            mutated_parameters.append((n, p, high_vert))
        if low_vert > 0.5:
            mutated_parameters.append((n, p, low_vert))
    return mutated_parameters


def assess_fit(masses, n, p, vert, neg_penalty):
    """
    Returns a lower-is-better value quantifying how well the empirical distribution matches the
    gamma distribution. A value of 0 means a perfect fit.
    """
    negbin_x = np.arange(len(masses))
    negbin_y = nbinom.pmf(negbin_x, n, p)
    negbin_y = negbin_y * vert

    differences = masses - negbin_y

    # Positive differences are where the empirical distribution is higher than the gamma.
    pos_differences = np.clip(differences, 0.0, None)

    # Negative differences are where the empirical distribution is lower than the gamma.
    neg_differences = np.clip(differences, None, 0.0)

    return abs(neg_penalty * np.sum(neg_differences)) + np.sum(pos_differences)

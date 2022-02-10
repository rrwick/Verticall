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
from scipy.stats import gamma


def get_mean(masses):
    return np.average(range(len(masses)), weights=masses)


def get_variance(masses):
    mean = get_mean(masses)
    values = np.arange(len(masses))
    return np.average((values - mean)**2, weights=masses)


def get_shape_and_scale(mean, variance):
    """
    Converts the intuitive mean and variance metrics to the less-intuitive shape and scale
    parameters needed to define a gamma distribution.
    """
    shape = (mean * mean) / variance
    scale = variance / mean
    return shape, scale


def fit_gamma_to_distribution(masses):
    """
    This function takes an empirical discrete distribution and tries to find a gamma distribution
    which best fits the data.

    In addition to shape and scale (the gamma parameters), this function also adjusts the vertical
    scale (from 0.5 to 1.0) because the empirical distribution may be multimodal, and we want to
    fit the gamma distribution to only the strongest mode, which we assume is at least half of the
    distribution.

    Empirical points below the gamma are penalised more than empirical points above the gamma. This
    is because in a multimodal distribution we expect other peaks to rise above the strongest mode.
    """
    masses = np.array(masses)
    neg_penalty = 2.0  # Controls how much being low is worse than being high.

    # We start with shape/scale based on the mean/variance of the distribution and a vertical scale
    # of one. For unimodal distributions, this should get use pretty close to the correct answer.
    mean, variance = get_mean(masses), get_variance(masses)
    shape, scale = get_shape_and_scale(mean, variance)
    vert = 1.0
    best_shape, best_scale, best_vert = shape, scale, vert
    best_fit = assess_fit(masses, shape, scale, vert, neg_penalty=neg_penalty)

    shape_step, scale_step, vert_step = 0.1, 0.1, 0.1
    while True:
        improvement = False
        mutated_parameters = get_mutated_parameters(best_shape, best_scale, best_vert,
                                                    shape_step, scale_step, vert_step)
        for shape, scale, vert in mutated_parameters:
            fit = assess_fit(masses, shape, scale, vert, neg_penalty=neg_penalty)
            if best_fit - fit > 0.0000001:
                best_shape, best_scale, best_vert = shape, scale, vert
                best_fit = fit
                improvement = True
        if not improvement:
            shape_step /= 2.0
            scale_step /= 2.0
            vert_step /= 2.0
        if shape_step < 0.0000001:
            break

    return best_shape, best_scale, best_vert


def get_mutated_parameters(shape, scale, vert, shape_step, scale_step, vert_step):
    mutated_parameters = []
    if shape_step > 0.0:
        high_shape = shape + shape_step
        low_shape = shape - shape_step
        mutated_parameters.append((high_shape, scale, vert))
        if low_shape > 0.0:
            mutated_parameters.append((low_shape, scale, vert))
    if scale_step > 0.0:
        high_scale = scale + scale_step
        low_scale = scale - scale_step
        mutated_parameters.append((shape, high_scale, vert))
        if low_scale > 0.0:
            mutated_parameters.append((shape, low_scale, vert))
    if vert_step > 0.0:
        high_vert = vert + vert_step
        low_vert = vert - vert_step
        if high_vert <= 1.0:
            mutated_parameters.append((shape, scale, high_vert))
        if low_vert > 0.5:
            mutated_parameters.append((shape, scale, low_vert))
    return mutated_parameters


def assess_fit(masses, shape, scale, vert, neg_penalty):
    """
    Returns a lower-is-better value quantifying how well the empirical distribution matches the
    gamma distribution. A value of 0 means a perfect fit.
    """
    gamma_x = np.arange(len(masses))
    gamma_y = gamma.pdf(gamma_x, shape, scale=scale)
    gamma_y = gamma_y * vert

    differences = masses - gamma_y

    # Positive differences are where the empirical distribution is higher than the gamma.
    pos_differences = np.clip(differences, 0.0, None)

    # Negative differences are where the empirical distribution is lower than the gamma.
    neg_differences = np.clip(differences, None, 0.0)

    return abs(neg_penalty * np.sum(neg_differences)) + np.sum(pos_differences)

# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 Davide Fioriti
#
# SPDX-License-Identifier: GPL-3.0-or-later
# coding: utf-8

module PointEstimateMethod

using Distributions
using Combinatorics
using JuMP
using HiGHS
using Polynomials

export pem

DEFAULT_SOLVER = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

include("auxiliaries.jl")
include("pem.jl")

end

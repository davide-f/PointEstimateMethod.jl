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
import Polynomials

export pem

DEFAULT_SOLVER = optimizer_with_attributes(
    HiGHS.Optimizer,
    "output_flag" => false,
    "primal_feasibility_tolerance"=>1e-06,
    "dual_feasibility_tolerance"=>1e-06,
    "ipm_optimality_tolerance"=>1e-07,
)

include("auxiliaries.jl")
include("pem.jl")

end

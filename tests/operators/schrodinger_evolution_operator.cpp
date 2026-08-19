/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include "catch2/catch_all.hpp"

#include "factory_functions.h"

#include "functions/GaussFunc.h"
#include "functions/special_functions.h"
#include "operators/MWOperator.h"
#include "operators/TimeEvolutionOperator.h"
#include "treebuilders/add.h"
#include "treebuilders/apply.h"
#include "treebuilders/project.h"

namespace schrodinger_evolution_operator {

TEST_CASE("Apply Schrodinger's evolution operator", "[apply_schrodinger_evolution], [schrodinger_evolution_operator], [mw_operator]") {
    const auto min_scale = 0;
    const auto max_depth = 25;

    const auto order = 4;
    const auto prec = 1.0e-7;

    int finest_scale = 7; // for time evolution operator construction (not recommended to use more than 10)

    // Time moments:
    double t1 = 0.001;        // initial time moment (not recommended to use more than 0.001)
    double delta_t = 0.03;    // time step (not recommended to use less than 0.001)
    double t2 = delta_t + t1; // final time moment

    // Initialize world in the unit cube [0,1]
    auto basis = mrcpp::LegendreBasis(order);
    auto world = mrcpp::BoundingBox<1>(min_scale);
    auto MRA = mrcpp::MultiResolutionAnalysis<1>(world, basis, max_depth);

    // Time evolution operator Exp(delta_t): one complex kernel, no Re/Im split
    mrcpp::TimeEvolutionOperator<1> Exp(MRA, prec, delta_t, finest_scale);

    // Analytical solution parameters for psi(x, t)
    double sigma = 0.001;
    double x0 = 0.5;

    // Functions f(x) = psi(x, t1) and g(x) = psi(x, t2)
    auto f = [sigma, x0, t = t1](const mrcpp::Coord<1> &r) -> ComplexDouble { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma); };
    auto g = [sigma, x0, t = t2](const mrcpp::Coord<1> &r) -> ComplexDouble { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma); };

    // Projecting functions
    mrcpp::FunctionTree<1, ComplexDouble> f_tree(MRA);
    mrcpp::project<1, ComplexDouble>(prec, f_tree, f);
    mrcpp::FunctionTree<1, ComplexDouble> g_tree(MRA);
    mrcpp::project<1, ComplexDouble>(prec, g_tree, g);

    // Apply operator Exp(delta_t) f(x)
    mrcpp::FunctionTree<1, ComplexDouble> fout_tree(MRA);
    mrcpp::apply<1, ComplexDouble>(prec, fout_tree, Exp, f_tree, -1, false);

    // Check g(x) = Exp(delta_t) f(x)
    mrcpp::FunctionTree<1, ComplexDouble> error(MRA); // = fout_tree - g_tree
    mrcpp::add<1, ComplexDouble>(prec, error, {1.0, 0.0}, fout_tree, {-1.0, 0.0}, g_tree, -1, false, false);
    auto sq_norm = error.getSquareNorm(); // Re: 1.7e-16, Im: 1.7e-17 when split

    double tolerance = prec * prec / 25.0; // 4.0e-16 (both components in one norm)

    REQUIRE(sq_norm == Catch::Approx(0.0).margin(tolerance));
}

} // namespace schrodinger_evolution_operator

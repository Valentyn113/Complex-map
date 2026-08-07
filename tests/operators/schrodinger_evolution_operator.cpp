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
#include "treebuilders/complex_apply.h"
#include "treebuilders/project.h"
#include "trees/MWNode.h"
#include "trees/OperatorTree.h"
#include "utils/math_utils.h"
 
#include <algorithm>
#include <cmath>
#include <complex>
 
using namespace mrcpp;
 
namespace schrodinger_evolution_operator {
 
/** @brief Parameters shared by the real-split and native-complex scenarios.
 *
 *  Keeping these in one place means the two constructions below are compared
 *  under identical conditions: same basis, same world, same time step and
 *  same finest scale for the operator build.
 */
namespace {
const auto min_scale = 0;
const auto max_depth = 25;
const auto order = 4;
const auto prec = 1.0e-7;
 
const int finest_scale = 7;  // for time evolution operator construction (not recommended to use more than 10)
const double delta_t = 0.03; // time step (not recommended to use less than 0.001)
 
// Analytical free-particle solution parameters
const double sigma = 0.001;
const double x0 = 0.5;
const double t1 = 0.001;        // initial time moment (not recommended to use more than 0.001)
const double t2 = delta_t + t1; // final time moment
 
auto make_MRA() {
    auto basis = LegendreBasis(order);
    auto world = BoundingBox<1>(min_scale);
    return MultiResolutionAnalysis<1>(world, basis, max_depth);
}
} // namespace
 
/** @brief Apply the Schrodinger semigroup to a free particle and compare with the
 *  analytical solution.
 *
 *  This is the classic construction: exp(i t d^2) is represented as two *real*
 *  operator trees, imaginary = false giving the cos part and imaginary = true the
 *  sin part, applied through the 2x2 real block via ComplexObject. This path is
 *  what all current callers use, so it is kept under test in its own right.
 */
TEST_CASE("Apply Schrodinger's evolution operator", "[apply_schrodinger_evolution], [schrodinger_evolution_operator], [mw_operator]") {
    auto MRA = make_MRA();
 
    // Time evolution operatror Exp(delta_t)
    mrcpp::TimeEvolutionOperator<1> ReExp(MRA, prec, delta_t, finest_scale, false);
    mrcpp::TimeEvolutionOperator<1> ImExp(MRA, prec, delta_t, finest_scale, true);
 
    // Functions f(x) = psi(x, t1) and g(x) = psi(x, t2)
    auto Re_f = [t = t1](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).real(); };
    auto Im_f = [t = t1](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).imag(); };
    auto Re_g = [t = t2](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).real(); };
    auto Im_g = [t = t2](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).imag(); };
 
    // Projecting functions
    mrcpp::FunctionTree<1> Re_f_tree(MRA);
    mrcpp::project<1, double>(prec, Re_f_tree, Re_f);
    mrcpp::FunctionTree<1> Im_f_tree(MRA);
    mrcpp::project<1, double>(prec, Im_f_tree, Im_f);
    mrcpp::FunctionTree<1> Re_g_tree(MRA);
    mrcpp::project<1, double>(prec, Re_g_tree, Re_g);
    mrcpp::FunctionTree<1> Im_g_tree(MRA);
    mrcpp::project<1, double>(prec, Im_g_tree, Im_g);
 
    // Output function trees
    mrcpp::FunctionTree<1> Re_fout_tree(MRA);
    mrcpp::FunctionTree<1> Im_fout_tree(MRA);
 
    // Complex objects for use in apply()
    mrcpp::ComplexObject<mrcpp::ConvolutionOperator<1>> E(ReExp, ImExp);
    mrcpp::ComplexObject<mrcpp::FunctionTree<1>> input(Re_f_tree, Im_f_tree);
    mrcpp::ComplexObject<mrcpp::FunctionTree<1>> output(Re_fout_tree, Im_fout_tree);
 
    // Apply operator Exp(delta_t) f(x)
    mrcpp::apply(prec, output, E, input);
 
    // Check g(x) = Exp(delta_t) f(x)
    mrcpp::FunctionTree<1> Re_error(MRA); // = Re_fout_tree - Re_g_tree
    mrcpp::FunctionTree<1> Im_error(MRA); // = Im_fout_tree - Im_g_tree
 
    // Re_error = Re_fout_tree - Re_g_tree
    mrcpp::add(prec, Re_error, 1.0, Re_fout_tree, -1.0, Re_g_tree);
    auto Re_sq_norm = Re_error.getSquareNorm(); // 1.7e-16
 
    // Im_error = Im_fout_tree - Im_g_tree
    mrcpp::add(prec, Im_error, 1.0, Im_fout_tree, -1.0, Im_g_tree);
    auto Im_sq_norm = Im_error.getSquareNorm(); // 1.7e-17
 
    double tolerance = prec * prec / 50.0; // 2.0e-16
 
    REQUIRE(Re_sq_norm == Catch::Approx(0.0).margin(tolerance));
    REQUIRE(Im_sq_norm == Catch::Approx(0.0).margin(tolerance));
}
 
/** @brief The native complex Schrodinger kernel must equal the real cos + i*sin split.
 *
 *  The classic construction above builds the Schrodinger semigroup exp(i t d^2) as
 *  two *real* operator trees. TimeEvolutionOperator<1, ComplexDouble> instead stores
 *  a single *complex* kernel whose coefficients are cos + i*sin.
 *
 *  Because the uniform tree builder fixes the operator-tree structure from the
 *  scale alone (independently of the coefficient values), the three operator trees
 *  are structurally identical, so we can compare them leaf by leaf. Every step that
 *  produces the coefficients (the cross-correlation contraction, the per-node MW
 *  compression, and the bottom-up MW transform) is linear with real filters, so the
 *  complex coefficients must equal (real coef) + i*(imag coef) to round-off.
 */
TEST_CASE("Complex Schrodinger evolution operator equals the real/imag split",
          "[complex_schrodinger_evolution], [schrodinger_evolution_operator], [mw_operator]") {
    auto MRA = make_MRA();
 
    // Classic split: two real operator trees (cos and sin parts).
    TimeEvolutionOperator<1, double> ReExp(MRA, prec, delta_t, finest_scale, false);
    TimeEvolutionOperator<1, double> ImExp(MRA, prec, delta_t, finest_scale, true);
 
    // Native single complex kernel (the bool flag is ignored for the complex type).
    TimeEvolutionOperator<1, ComplexDouble> Exp(MRA, prec, delta_t, finest_scale, false);
 
    const OperatorTree<double> &re = ReExp.getComponent(0, 0);
    const OperatorTree<double> &im = ImExp.getComponent(0, 0);
    const OperatorTree<ComplexDouble> &cx = Exp.getComponent(0, 0);
 
    // (1) Square-norm identity: ||cos + i sin||^2 == ||cos||^2 + ||sin||^2.
    const double sq_re = re.getSquareNorm();
    const double sq_im = im.getSquareNorm();
    const double sq_cx = cx.getSquareNorm();
    REQUIRE(sq_cx == Catch::Approx(sq_re + sq_im).epsilon(1.0e-10));
 
    // (2) Leaf-by-leaf: complex coef == real coef + i * imag coef.
    REQUIRE(cx.getNEndNodes() == re.getNEndNodes());
    REQUIRE(cx.getNEndNodes() == im.getNEndNodes());
 
    double max_diff = 0.0;
    for (int n = 0; n < cx.getNEndNodes(); n++) {
        const MWNode<2, ComplexDouble> &cx_node = cx.getEndMWNode(n);
        const MWNode<2, double> &re_node = re.getEndMWNode(n);
        const MWNode<2, double> &im_node = im.getEndMWNode(n);
 
        // The uniform builder must produce identical structure in all three trees.
        REQUIRE(cx_node.getNodeIndex() == re_node.getNodeIndex());
        REQUIRE(cx_node.getNodeIndex() == im_node.getNodeIndex());
 
        const int n_coefs = cx_node.getNCoefs();
        REQUIRE(re_node.getNCoefs() == n_coefs);
        REQUIRE(im_node.getNCoefs() == n_coefs);
 
        const ComplexDouble *cx_c = cx_node.getCoefs();
        const double *re_c = re_node.getCoefs();
        const double *im_c = im_node.getCoefs();
        for (int k = 0; k < n_coefs; k++) {
            const ComplexDouble expected(re_c[k], im_c[k]);
            max_diff = std::max(max_diff, std::abs(cx_c[k] - expected));
        }
    }
 
    REQUIRE(max_diff == Catch::Approx(0.0).margin(1.0e-12));
}
 
/** @brief End-to-end: evolve a free particle with the single complex operator.
 *
 *  The first scenario evolves psi by splitting exp(i t d^2) into two real operators
 *  and recombining them through the 2x2 ComplexObject block. Here the same evolution
 *  is done in one shot: a complex wavefunction in a single FunctionTree<1, ComplexDouble>,
 *  acted on by a single TimeEvolutionOperator<1, ComplexDouble>, compared against the
 *  analytical solution at the later time.
 *
 *  This exercises the whole complex path -- operator construction, the convolution with
 *  a complex operator, and apply -- rather than only the operator coefficients.
 */
TEST_CASE("Apply the complex Schrodinger evolution operator",
          "[apply_complex_schrodinger_evolution], [schrodinger_evolution_operator], [mw_operator]") {
    auto MRA = make_MRA();
 
    // Single complex kernel: exp(i * delta_t * d^2/dx^2)
    TimeEvolutionOperator<1, ComplexDouble> Exp(MRA, prec, delta_t, finest_scale, false);
 
    // psi(x, t1) and the reference psi(x, t2), both as genuine complex functions
    auto psi_t1 = [t = t1](const Coord<1> &r) -> ComplexDouble { return free_particle_analytical_solution(r[0], x0, t, sigma); };
    auto psi_t2 = [t = t2](const Coord<1> &r) -> ComplexDouble { return free_particle_analytical_solution(r[0], x0, t, sigma); };
 
    FunctionTree<1, ComplexDouble> f_tree(MRA);
    project<1, ComplexDouble>(prec, f_tree, psi_t1);
 
    FunctionTree<1, ComplexDouble> ref_tree(MRA);
    project<1, ComplexDouble>(prec, ref_tree, psi_t2);
 
    // psi(t2) = Exp(delta_t) psi(t1), in a single application
    FunctionTree<1, ComplexDouble> out_tree(MRA);
    apply(prec, out_tree, Exp, f_tree);
 
    // error = out - ref
    FunctionTree<1, ComplexDouble> error(MRA);
    add(prec, error, ComplexDouble(1.0, 0.0), out_tree, ComplexDouble(-1.0, 0.0), ref_tree);
 
    // getSquareNorm() on a complex tree sums |c|^2, i.e. both parts at once, so allow
    // twice the per-part tolerance used in the real/imag split scenario above.
    const double tolerance = 2.0 * prec * prec / 50.0;
    REQUIRE(error.getSquareNorm() == Catch::Approx(0.0).margin(tolerance));
}
 
} // namespace schrodinger_evolution_operator
 
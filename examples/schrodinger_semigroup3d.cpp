#include "MRCPP/MWFunctions"
#include "functions/special_functions.h"
#include "operators/TimeEvolutionOperator.h"
#include "treebuilders/apply.h"
#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include <cmath>
#include <cstdio>
#include <memory>

const auto min_scale = 0;
const auto max_depth = 25;

const auto order = 4;
const auto prec = 1.0e-4;

const auto D = 3;

// Time moments:
double t1 = 0.0005;      // initial time moment
double delta_t = 0.0005; // time step
int n_steps = 4;         // number of steps, the final time moment is t1 + n_steps * delta_t

/**
 * @brief Free-particle time evolution in 3D.
 * @details The operator is built once and applied at every step:
 * \f[
 *   \psi(\mathbf{r}, t + \delta t) = \exp \left( i \delta t \nabla^2 \right) \psi(\mathbf{r}, t)
 *   .
 * \f]
 *
 * The 3D propagator is a product of 1D propagators, one in each direction,
 * and the operator uses the same tree for x, y and z. A product of 1D
 * Gaussians therefore stays a product, and every step is compared with
 * \f[
 *   \psi(\mathbf{r}, t) = \prod_{d=1}^{3} \sqrt{\frac{\sigma}{4it + \sigma}} e^{-\frac{(x_d - x_0)^2}{4it + \sigma}}
 *   .
 * \f]
 *
 * The operator is the free-particle propagator restricted to the unit cube,
 * so the comparison holds as long as psi stays negligible at the boundary.
 *
 */
int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    // Initialize printing
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::print::environment(0);

    // Initialize world in the unit cube [0,1]^3
    auto basis = mrcpp::LegendreBasis(order);
    auto world = mrcpp::BoundingBox<D>(min_scale);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);

    mrcpp::print::header(0, "Building operator");

    // Time evolution operator Exp(delta_t)
    mrcpp::TimeEvolutionOperator<D> Exp(MRA, prec, delta_t);

    mrcpp::print::footer(0, timer, 2);

    // Analytical solution parameters for psi(r, t)
    double sigma = 0.01;
    double x0 = 0.5;
    auto psi = [sigma, x0](const mrcpp::Coord<D> &r, double t) -> ComplexDouble {
        ComplexDouble v = 1.0;
        for (int d = 0; d < D; d++) v *= mrcpp::free_particle_analytical_solution(r[d], x0, t, sigma);
        return v;
    };

    mrcpp::print::header(0, "Propagating");

    // Initial function psi(r, t1)
    auto f_tree = std::make_unique<mrcpp::FunctionTree<D, ComplexDouble>>(MRA);
    mrcpp::project<D, ComplexDouble>(prec, *f_tree, [&](const mrcpp::Coord<D> &r) { return psi(r, t1); });

    for (int n = 1; n <= n_steps; n++) {
        double t = t1 + n * delta_t;

        // One step: psi(r, t) = Exp(delta_t) psi(r, t - delta_t)
        auto step_timer = mrcpp::Timer();
        auto g_tree = std::make_unique<mrcpp::FunctionTree<D, ComplexDouble>>(MRA);
        mrcpp::apply<D, ComplexDouble>(prec, *g_tree, Exp, *f_tree, -1, false);
        step_timer.stop();
        f_tree = std::move(g_tree);

        // Compare with the analytical solution at time t
        mrcpp::FunctionTree<D, ComplexDouble> exact(MRA);
        mrcpp::project<D, ComplexDouble>(prec, exact, [&](const mrcpp::Coord<D> &r) { return psi(r, t); });
        mrcpp::FunctionTree<D, ComplexDouble> error(MRA);
        mrcpp::add<D, ComplexDouble>(prec, error, {1.0, 0.0}, *f_tree, {-1.0, 0.0}, exact, -1, false, false);

        char label[32];
        std::snprintf(label, sizeof(label), "Step %d, t = %.4f", n, t);
        mrcpp::print::tree(0, label, f_tree->getNNodes(), f_tree->getSizeNodes(), step_timer.elapsed());
        mrcpp::print::value(0, "Relative error", std::sqrt(error.getSquareNorm() / exact.getSquareNorm()));
    }

    mrcpp::print::footer(0, timer, 2);
    return 0;
}

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

#pragma once

#include "ConvolutionOperator.h"
#include "MWOperator.h"
#include "core/SchrodingerEvolution_CrossCorrelation.h"

namespace mrcpp {

/** @class TimeEvolutionOperator
 *
 * @brief Semigroup of the free-particle Schrodinger equation
 *
 * @details Represents the semigroup
 * \f$
 *      \exp \left( i t \partial_x^2 \right)
 *      .
 * \f$
 * Matrix elements (actual operator tree) of the operator can be obtained by calling getComponent(0, 0).
 *
 * @note So far implementation is done for Legendre scaling functions in 1d.
 *
 * \todo: Extend to D dimensinal on a general interval [a, b] in the future.
 *
 */
template <int D>
class TimeEvolutionOperator : public ConvolutionOperator<D> // One can use ConvolutionOperator instead as well
{
public:
    /** @brief Semigroup \f$ \exp(i t \partial_x^2) \f$ at one time moment.
     *
     * @param[in] mra: which MRA to operate on
     * @param[in] prec: build precision
     * @param[in] time: time step \f$ t \f$
     *
     * @details Built adaptively from the precision, as every other convolution
     * operator is. The kernel is complex and is held as a single
     * `OperatorTree<ComplexDouble>` -- the operator-side counterpart of a
     * complex `FunctionTree` -- with refinement thresholding on the modulus of
     * the coefficients.
     *
     * @note Applying this to a real `FunctionTree` aborts; project the input as
     * `ComplexDouble` first, or let the `CompFunction` overload promote it.
     */
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time);

    /** @brief As above, refined uniformly down to a given scale.
     *
     * @param[in] finest_scale: uniform refinement down to this scale
     * @param[in] max_Jpower: number of power integrals retained
     *
     * @details Bypasses the adaptive build. Kept for callers that need a fixed
     * grid, in the same spirit as the `root`/`reach` overloads of
     * `PoissonOperator` and `HelmholtzOperator`.
     */
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, int max_Jpower = 30);

    /** @brief Rejects the old `bool imaginary` argument at compile time.
     *
     * @details The semigroup is no longer built one part at a time. A `bool` in
     * that position would otherwise bind silently to `finest_scale`.
     */
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, bool imaginary, int max_Jpower = 30) = delete;

    /** @brief Rejects the old fixed-scale signature at compile time.
     *
     * @details `(mra, prec, time, 7, false)` would otherwise bind to the
     * `finest_scale, max_Jpower` overload with `max_Jpower = 0`, since `int`
     * is an exact match in the fourth position and `bool` converts to `int`
     * in the fifth.
     */
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, bool imaginary, int max_Jpower = 30) = delete;
    TimeEvolutionOperator(const TimeEvolutionOperator &oper) = delete;
    TimeEvolutionOperator &operator=(const TimeEvolutionOperator &oper) = delete;
    virtual ~TimeEvolutionOperator() = default;

    double getBuildPrec() const { return this->build_prec; }

protected:
    void initialize(double time, int finest_scale, int max_Jpower);
    void initialize(double time, int max_Jpower);
    void initializeSemiUniformly(double time, int max_Jpower);

    void setBuildPrec(double prec) { this->build_prec = prec; }

    double build_prec{-1.0};
    SchrodingerEvolution_CrossCorrelation *cross_correlation{nullptr};
};

} // namespace mrcpp

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

#include "TreeCalculator.h"
#include "core/CrossCorrelationCache.h"
#include "core/SchrodingerEvolution_CrossCorrelation.h"
#include "functions/JpowerIntegrals.h"

namespace mrcpp {

/** @class TimeEvolution_CrossCorrelationCalculator
 *
 * @brief An efficient way to calculate ... (work in progress)
 *
 * @details An efficient way to calculate ... having the form
 * \f$ \ldots = \ldots \f$
 *
 *
 *
 */
/** @brief Which part of the complex kernel a calculator pass writes.
 *
 * @details `Kernel_Modulus` writes \f$ |z| \f$ and exists only to drive
 * adaptive refinement: the grid is settled once on the modulus, then the two
 * real parts are filled on it. This is the operator-side counterpart of
 * `FunctionTree::Real()` and `FunctionTree::Imag()`, which likewise derive
 * real trees from the grid of the complex object rather than refining each
 * part on its own.
 */
enum Kernel_Part { Kernel_Real, Kernel_Imag, Kernel_Modulus };

class TimeEvolution_CrossCorrelationCalculator final : public TreeCalculator<2> {
public:
    TimeEvolution_CrossCorrelationCalculator(std::map<int, JpowerIntegrals *> &J, SchrodingerEvolution_CrossCorrelation *cross_correlation, Kernel_Part part)
            : J_power_inetgarls(J)
            , cross_correlation(cross_correlation)
            , part(part) {}
    // private:
    std::map<int, JpowerIntegrals *> J_power_inetgarls;
    SchrodingerEvolution_CrossCorrelation *cross_correlation;

    /// @brief Part of the complex kernel written by this pass.
    Kernel_Part part;

    void calcNode(MWNode<2> &node) override;

    // template <int T>
    void applyCcc(MWNode<2> &node);
    // template <int T> void applyCcc(MWNode<2> &node, CrossCorrelationCache<T> &ccc);
};

} // namespace mrcpp

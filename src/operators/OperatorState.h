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

/** OperatorState is a simple helper class for operator application.
 * It keeps track of various state dependent variables and memory
 * regions. We cannot have some of this information directly in OperatorFunc
 * because of multi-threading issues.
 */

#pragma once

#include <vector>

#include <Eigen/Core>

#include "trees/MWNode.h"
#include "utils/math_utils.h"

namespace mrcpp {

namespace detail {

/** @brief Map one operator band with the scalar type of the function tree it multiplies.
 *
 * Operator coefficients are stored as ComplexDouble unconditionally. A real
 * function tree can only be convolved with a real kernel, so for T = double we
 * map the real parts in place with an inner stride of 2; this is well defined
 * because std::complex<double> is layout-compatible with double[2] ([complex.numbers.general]).
 * Flop count on the real path is therefore unchanged; only the load pattern is strided.
 */
template <typename T> struct OperBandMap;

template <> struct OperBandMap<ComplexDouble> {
    using Type = Eigen::Map<const Eigen::MatrixXcd>;
    static Type get(const ComplexDouble *o, int kp1) { return Type(o, kp1, kp1); }
};

template <> struct OperBandMap<double> {
    using StrideType = Eigen::Stride<Eigen::Dynamic, 2>;
    using Type = Eigen::Map<const Eigen::MatrixXd, Eigen::Unaligned, StrideType>;
    static Type get(const ComplexDouble *o, int kp1) {
        return Type(reinterpret_cast<const double *>(o), kp1, kp1, StrideType(2 * kp1, 2));
    }
};

} // namespace detail

#define GET_OP_IDX(FT, GT, ID) (2 * ((GT >> ID) & 1) + ((FT >> ID) & 1))

template <int D, typename T> class OperatorState final {
public:
    OperatorState(MWNode<D, T> &gn, T *scr1)
            : gNode(&gn) {
        this->kp1 = this->gNode->getKp1();
        this->kp1_d = this->gNode->getKp1_d();
        this->kp1_2 = math_utils::ipow(this->kp1, 2);
        this->kp1_dm1 = math_utils::ipow(this->kp1, D - 1);
        this->gData = this->gNode->getCoefs();
        this->maxDeltaL = -1;

        T *scr2 = scr1 + this->kp1_d;

        for (int i = 1; i < D; i++) {
            if (IS_ODD(i)) {
                this->aux[i] = scr2;
            } else {
                this->aux[i] = scr1;
            }
        }
    }

    OperatorState(MWNode<D, T> &gn, std::vector<T> scr1)
            : OperatorState(gn, scr1.data()) {}
    void setFNode(MWNode<D, T> &fn) {
        this->fNode = &fn;
        this->fData = this->fNode->getCoefs();
    }
    void setFIndex(NodeIndex<D> &idx) {
        this->fIdx = &idx;
        calcMaxDeltaL();
    }
    void setGComponent(int gt) {
        this->aux[D] = this->gData + gt * this->kp1_d;
        this->gt = gt;
    }
    void setFComponent(int ft) {
        this->aux[0] = this->fData + ft * this->kp1_d;
        this->ft = ft;
    }

    int getMaxDeltaL() const { return this->maxDeltaL; }
    int getOperIndex(int i) const { return GET_OP_IDX(this->ft, this->gt, i); }

    T **getAuxData() { return this->aux; }
    ComplexDouble **getOperData() { return this->oData; }

    friend class ConvolutionCalculator<D, T>;
    friend class DerivativeCalculator<D, T>;

private:
    int ft;
    int gt;

    int maxDeltaL;
    double fThreshold;
    double gThreshold;
    // Shorthands
    int kp1;
    int kp1_2;
    int kp1_d;
    int kp1_dm1;

    MWNode<D, T> *gNode;
    MWNode<D, T> *fNode;
    NodeIndex<D> *fIdx;

    T *aux[D + 1];
    T *gData;
    T *fData;
    ComplexDouble *oData[D];

    void calcMaxDeltaL() {
        const auto &gl = this->gNode->getNodeIndex();
        const auto &fl = *this->fIdx;
        int max_dl = 0;
        for (int d = 0; d < D; d++) {
            int dl = std::abs(fl[d] - gl[d]);
            if (dl > max_dl) { max_dl = dl; }
        }
        this->maxDeltaL = max_dl;
    }
};

} // namespace mrcpp

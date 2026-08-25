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

#include "TimeEvolutionOperator.h"
//#include "MRCPP/MWOperators"

#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"

#include "functions/GaussExp.h"
#include "functions/Gaussian.h"

#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/DefaultCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/SplitAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"

#include "trees/BandWidth.h"
#include "trees/CornerOperatorTree.h"
#include "trees/FunctionTreeVector.h"

#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"

#include "treebuilders/TimeEvolution_CrossCorrelationCalculator.h"

#include <vector>

#include "trees/OperatorNode.h"

namespace mrcpp {

namespace {

/** @brief True where an already-built reference tree has children. */
inline bool refines_in(const OperatorTree &reference, const MWNode<2> &node) {
    const MWNode<2> *ref_node = const_cast<OperatorTree &>(reference).findNode(node.getNodeIndex());
    return (ref_node != nullptr) and ref_node->isBranchNode();
}

/** @brief The splitting rule of `OperatorAdaptor`.
 *
 * @details Duplicated because `OperatorAdaptor` is `final` and cannot be
 * derived from.
 */
inline bool operator_split_check(const MWNode<2> &node) {
    const auto &idx = node.getNodeIndex();
    bool chkTransl = (idx[0] == 0 or idx[1] == 0);
    bool chkCompNorm = false;
    for (int i = 1; i < 4; i++) {
        if (node.getComponentNorm(i) > 0.0) chkCompNorm = true;
    }
    return chkTransl and chkCompNorm;
}

/** @class GridCopyAdaptor
 *
 * @brief Refines a tree onto the grid of an already-built reference tree.
 *
 * @details `CopyAdaptor` does the same job for `FunctionTree`, but takes a
 * `FunctionTreeVector` and so cannot be pointed at an `OperatorTree`.
 */
class GridCopyAdaptor final : public TreeAdaptor<2> {
public:
    GridCopyAdaptor(const OperatorTree &ref, int ms)
            : TreeAdaptor<2>(ms)
            , reference(ref) {}

protected:
    const OperatorTree &reference;

    bool splitNode(const MWNode<2> &node) const override { return refines_in(this->reference, node); }
};

/** @class GridUnionAdaptor
 *
 * @brief Splits where the operator criterion fires or the reference tree does.
 *
 * @details Used for the second of three build passes, so that the resulting
 * grid is the union of what the two halves of the kernel need.
 */
class GridUnionAdaptor final : public TreeAdaptor<2> {
public:
    GridUnionAdaptor(const OperatorTree &ref, int ms)
            : TreeAdaptor<2>(ms)
            , reference(ref) {}

protected:
    const OperatorTree &reference;

    bool splitNode(const MWNode<2> &node) const override { return operator_split_check(node) or refines_in(this->reference, node); }
};

} // namespace


/** @brief A uniform constructor for TimeEvolutionOperator class.
 *
 * @param[in] mra: MRA.
 * @param[in] prec: precision.
 * @param[in] time: the time moment (step).
 * @param[in] finest_scale: the operator tree is constructed uniformly down to this scale.
 * @param[in] max_Jpower: maximum amount of power integrals used.
 *
 * @details Constructs the Schrodinger semigroup at a given time moment, as a matched pair of real operator trees holding its real and imaginary parts.
 *
 */
template <int D>
TimeEvolutionOperator<D>::TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, int max_Jpower)
        : ConvolutionOperator<D>(mra, mra.getRootScale(), -10) // One can use ConvolutionOperator instead as well
{
    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);

    SchrodingerEvolution_CrossCorrelation cross_correlation(30, mra.getOrder(), mra.getScalingBasis().getScalingType());
    this->cross_correlation = &cross_correlation;

    // will go outside of the constructor in future
    if (finest_scale < 0)
        initialize(time, max_Jpower);
    else
        initialize(time, finest_scale, max_Jpower);

    this->initOperExp(1); // this turns out to be important
    Printer::setPrintLevel(oldlevel);
}

/** @brief An adaptive constructor for TimeEvolutionOperator class.
 *
 * @param[in] mra: MRA.
 * @param[in] prec: precision.
 * @param[in] time: the time moment (step).
 * @param[in] max_Jpower: maximum amount of power integrals used.
 *
 * @details Adaptively constructs the Schrodinger semigroup at a given time moment, as a matched pair of real operator trees.
 * It is recommended for use in case of high polynomial order in use of the scaling basis.
 *
 * @note For technical reasons the operator tree is constructed no deeper than to scale \f$ n = 18 \f$.
 * This should be weakened in future.
 *
 */

/** @brief Creates the operator
 *
 * @details Adaptive down to scale \f$ N = 18 \f$.
 * This scale limit bounds the amount of JpowerIntegrals
 * to be calculated.
 * @note In future work we plan to optimize calculation of JpowerIntegrals so that we calculate
 * only needed ones, while building the tree (in progress).
 *
 */
template <int D> void TimeEvolutionOperator<D>::initialize(double time, int max_Jpower) {
    int N = 18;

    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();

    std::map<int, JpowerIntegrals *> J;
    for (int n = 0; n <= N + 1; n++) J[n] = new JpowerIntegrals(time * std::pow(4, n), n, max_Jpower);

    mrcpp::TreeBuilder<2> builder;
    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale(), true);

    // The two halves must share a grid, since the contraction indexes both
    // trees with the same node, and each half needs whatever its own wavelet
    // norms demand. Take the union in three passes:
    //  - refine the real part on its own criterion
    //  - refine the imaginary part on its criterion OR the real grid
    //  - rebuild the real part on that union
    auto build_part = [&](Kernel_Part part, TreeAdaptor<2> &adapt) {
        auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);
        TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, part);
        builder.build(*o_tree, calculator, adapt, N);

        Timer trans_t;
        o_tree->mwTransform(BottomUp);
        o_tree->removeRoughScaleNoise();
        o_tree->calcSquareNorm();
        o_tree->setupOperNodeCache();
        print::time(10, "Time transform", trans_t);
        print::separator(10, ' ');
        return o_tree;
    };

    auto re_first = build_part(Kernel_Real, adaptor);

    GridUnionAdaptor union_adaptor(*re_first, o_mra.getMaxScale());
    auto im_tree = build_part(Kernel_Imag, union_adaptor);

    GridCopyAdaptor copy_adaptor(*im_tree, o_mra.getMaxScale());
    auto re_tree = build_part(Kernel_Real, copy_adaptor);

    this->raw_exp.push_back(std::move(re_tree));
    this->raw_exp_im.push_back(std::move(im_tree));

    for (int n = 0; n <= N + 1; n++) delete J[n];
}

/** @brief Creates the operator
 *
 * @details Uniform down to finest scale.
 *
 */
template <int D> void TimeEvolutionOperator<D>::initialize(double time, int finest_scale, int max_Jpower) {
    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();

    // Setup uniform tree builder
    TreeBuilder<2> builder;
    SplitAdaptor<2> uniform(o_mra.getMaxScale(), true);

    int N = finest_scale;
    double threshold = o_prec / 1000.0;
    std::map<int, JpowerIntegrals *> J;
    for (int n = 0; n <= N + 1; n++) J[n] = new JpowerIntegrals(time * std::pow(4, n), n, max_Jpower, threshold);

    // Uniform refinement, so the real and imaginary trees share a grid by
    // construction and the contraction can index both with the same node.
    auto build_part = [&](Kernel_Part part) {
        TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, part);
        auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);
        builder.build(*o_tree, calculator, uniform, N); // Expand 1D kernel into 2D operator

        Timer trans_t;
        o_tree->mwTransform(BottomUp);
        o_tree->calcSquareNorm();
        o_tree->setupOperNodeCache();
        print::time(10, "Time transform", trans_t);
        print::separator(10, ' ');
        return o_tree;
    };

    this->raw_exp.push_back(build_part(Kernel_Real));
    this->raw_exp_im.push_back(build_part(Kernel_Imag));

    for (int n = 0; n <= N + 1; n++) delete J[n];
}

/** @brief Creates the operator (in progress)
 *
 * @details Tree construction starts uniformly and then continues adaptively down to scale \f$ N = 18 \f$.
 * This scale limit bounds the amount of JpowerIntegrals
 * to be calculated.
 * @note This method is not ready for use and should not be used (in progress).
 *
 */
template <int D> void TimeEvolutionOperator<D>::initializeSemiUniformly(double time, int max_Jpower) {
    MSG_ERROR("Not implemented yet method.");

    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();

    mrcpp::TreeBuilder<2> builder;
    mrcpp::SplitAdaptor<2> uniform(o_mra.getMaxScale(), true);

    int N = 18;

    auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);
    DefaultCalculator<2> intitial_calculator;
    for (auto n = 0; n < 4; n++) builder.build(*o_tree, intitial_calculator, uniform, 1);

    double threshold = o_prec / 1000.0;
    std::map<int, mrcpp::JpowerIntegrals *> J;
    for (int n = 0; n <= N + 1; n++) J[n] = new mrcpp::JpowerIntegrals(time * std::pow(4, n), n, max_Jpower, threshold);
    mrcpp::TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, Kernel_Real);

    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale());
    builder.build(*o_tree, calculator, adaptor, 13);

    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(mrcpp::BottomUp);
    o_tree->removeRoughScaleNoise();
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; n++) delete J[n];
}

template class TimeEvolutionOperator<1>;
template class TimeEvolutionOperator<2>;
template class TimeEvolutionOperator<3>;

} // namespace mrcpp

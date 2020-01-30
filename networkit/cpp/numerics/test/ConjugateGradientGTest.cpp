/*
 * ConjugateGradientGTest.cpp
 *
 *  Created on: 30.01.2020
 *      Author: Eugenio Angrman <angrimae@hu-berlin.de>
 */

// networkit-format

#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/numerics/ConjugateGradient.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>
#include <networkit/numerics/Preconditioner/IdentityPreconditioner.hpp>

namespace NetworKit {

class ConjugateGradientGTest : public testing::Test {};

TEST_F(ConjugateGradientGTest, testSmallGraph) {
    METISGraphReader reader;
    const auto G = ConnectedComponents::extractLargestConnectedComponent(
        reader.read("input/power.graph"), true);
    const auto n = G.numberOfNodes();
    Vector rhs(n), resultLamg(n), resultCG(n);
    for (index i = 0; i < n; ++i) {
        rhs[i] = 2.0 * Aux::Random::probability() - 1.0;
        rhs[i] *= rhs[i];
    }

    const auto L = CSRMatrix::laplacianMatrix(G);
    static constexpr double tol = 1e-9;
    Lamg<CSRMatrix> solver(tol);
    solver.setupConnected(L);
    solver.solve(rhs, resultLamg);

    ConjugateGradient<CSRMatrix, IdentityPreconditioner> cg(tol);
    cg.setupConnected(L);
    cg.solve(rhs, resultCG);

    for (index i = 0; i < resultLamg.getDimension(); ++i) {
        INFO(resultLamg[i], " ", resultCG[i]);
    }
}
} // namespace NetworKit

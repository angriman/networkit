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
        reader.read("input/celegans_metabolic.graph"), true);
    const auto n = G.numberOfNodes();
    Vector rhs(n, 0), resultLamg(n), resultCG(n);
    node u = G.randomNode();
    node v = G.randomNode();
    while (v == u)
        v = G.randomNode();

    rhs[u] = 1;
    rhs[v] = -1;
    const auto L = CSRMatrix::laplacianMatrix(G);
    static constexpr double tol = 1e-6;
    Lamg<CSRMatrix> solver(tol);
    solver.setupConnected(L);
    solver.solve(rhs, resultLamg);
    INFO("LAMG SOLVED");

    ConjugateGradient<CSRMatrix, IdentityPreconditioner> cg(tol);
    cg.setupConnected(L);
    cg.solve(rhs, resultCG);

    for (index i = 0; i < resultLamg.getDimension(); ++i) {
        INFO(resultLamg[i], " ", resultCG[i]);
    }
}
} // namespace NetworKit

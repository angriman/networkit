// networkit-format
/*
 * EigenvectorCentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <cmath>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/EigenvectorCentrality.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>

#include "CentralityTestFixture.hpp"

namespace NetworKit {

class EigenvectorCentralityGTest : public DirectedWeightedGTest {};

INSTANTIATE_TEST_SUITE_P(InstantiationName, EigenvectorCentralityGTest, directedWeightedParams());

TEST_F(EigenvectorCentralityGTest, testEigenvectorCentrality) {
    /* Graph:
     0    3   6
      \  / \ /
       2 -- 5
      /  \ / \
     1    4   7

     Edges in the upper row have weight 3,
     the edge in the middle row has weight 1.5,
     edges in the lower row have weight 2.
    */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    EigenvectorCentrality centrality(G);
    centrality.run();
    std::vector<double> cen = centrality.scores();

    // computed with Matlab
    const double tol = 1e-4;
    EXPECT_NEAR(0.2254, std::fabs(cen[0]), tol);
    EXPECT_NEAR(0.1503, std::fabs(cen[1]), tol);
    EXPECT_NEAR(0.5290, std::fabs(cen[2]), tol);
    EXPECT_NEAR(0.4508, std::fabs(cen[3]), tol);
    EXPECT_NEAR(0.3006, std::fabs(cen[4]), tol);
    EXPECT_NEAR(0.5290, std::fabs(cen[5]), tol);
    EXPECT_NEAR(0.2254, std::fabs(cen[6]), tol);
    EXPECT_NEAR(0.1503, std::fabs(cen[7]), tol);
}

TEST_F(EigenvectorCentralityGTest, testEigenvectorCentralityEmptyGraph) {
    Graph G(0);

    EigenvectorCentrality centrality(G);
    centrality.run();

    EXPECT_TRUE(centrality.scores().empty());
}

TEST_P(EigenvectorCentralityGTest, testEigenvectorCentralityNoEdges) {
    for (count n : {1, 2, 5}) {
        Graph G(n, isWeighted(), isDirected());

        EigenvectorCentrality centrality(G);
        centrality.run();
        const std::vector<double> cen = centrality.scores();

        ASSERT_EQ(n, cen.size());

        double squaredNorm = 0.0;
        for (node u = 0; u < n; ++u) {
            EXPECT_TRUE(std::isfinite(cen[u]));
            EXPECT_NEAR(1.0 / std::sqrt(static_cast<double>(n)), cen[u], 1e-12);
            squaredNorm += cen[u] * cen[u];
        }
        EXPECT_NEAR(1.0, squaredNorm, 1e-12);
    }
}

TEST_F(EigenvectorCentralityGTest, testEigenvectorCentralityIsolatedNode) {
    count n = 9;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    EigenvectorCentrality centrality(G);
    centrality.run();
    std::vector<double> cen = centrality.scores();

    const double tol = 1e-4;
    EXPECT_NEAR(0.2254, std::fabs(cen[0]), tol);
    EXPECT_NEAR(0.1503, std::fabs(cen[1]), tol);
    EXPECT_NEAR(0.5290, std::fabs(cen[2]), tol);
    EXPECT_NEAR(0.4508, std::fabs(cen[3]), tol);
    EXPECT_NEAR(0.3006, std::fabs(cen[4]), tol);
    EXPECT_NEAR(0.5290, std::fabs(cen[5]), tol);
    EXPECT_NEAR(0.2254, std::fabs(cen[6]), tol);
    EXPECT_NEAR(0.1503, std::fabs(cen[7]), tol);
    EXPECT_NEAR(0.0, cen[8], tol);
}

TEST_F(EigenvectorCentralityGTest, testEigenvectorCentralityDirectedAcyclic) {
    Graph G(2, false, true);
    G.addEdge(0, 1);

    EigenvectorCentrality centrality(G);
    centrality.run();
    const std::vector<double> cen = centrality.scores();

    ASSERT_EQ(2u, cen.size());
    EXPECT_TRUE(std::isfinite(cen[0]));
    EXPECT_TRUE(std::isfinite(cen[1]));
}

TEST_F(EigenvectorCentralityGTest, benchEigenvectorCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    EigenvectorCentrality cen(G);
    cen.run();
    std::vector<std::pair<node, double>> ranking = cen.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

} /* namespace NetworKit */

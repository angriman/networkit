// networkit-format
/*
 * BetweennessGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class BetweennessGTest : public testing::Test {};

TEST_F(BetweennessGTest, testBetweennessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.scores();
    std::vector<double> result = {0., 0., 15., 3., 3., 1.};

    const double tol = 1e-3;
    G.forNodes([&](node u) { EXPECT_NEAR(result[u], bc[u], tol); });

    Betweenness centralityNorm(G, true);
    centralityNorm.run();
    bc = centralityNorm.scores();
    const double pairs = (n - 1) * (n - 2);
    G.forNodes([&](node u) { EXPECT_NEAR(result[u] / pairs, bc[u], tol); });
}

TEST_F(BetweennessGTest, testBetweenness2Centrality) {
    /* Graph:
          0    3
          \  / \
          2    5
          /  \ /
          1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    const double tol = 1e-3;
    EXPECT_NEAR(0.0, bc[0], tol);
    EXPECT_NEAR(0.0, bc[1], tol);
    EXPECT_NEAR(15.0, bc[2], tol);
    EXPECT_NEAR(3.0, bc[3], tol);
    EXPECT_NEAR(3.0, bc[4], tol);
    EXPECT_NEAR(1.0, bc[5], tol);
}

TEST_F(BetweennessGTest, testBetweennessCompactScores) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    G.removeNode(4); // remove a node to test compactScores

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.compactScores();

    Graph G2(n - 1);

    G2.addEdge(0, 2);
    G2.addEdge(1, 2);
    G2.addEdge(2, 3);
    G2.addEdge(3, 4);

    Betweenness centrality2(G);
    centrality2.run();
    std::vector<double> bc2 = centrality2.scores();
    for (size_t i = 0; i < bc.size(); i++) {
        EXPECT_EQ(bc[i], bc2[i]);
    }
}

TEST_F(BetweennessGTest, testBetweennessCompactEdgeScores) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.indexEdges();

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    G.removeEdge(2, 4); // remove an edge to test compactEdgeScores

    Betweenness centrality(G, false, true);
    centrality.run();
    std::vector<double> bc = centrality.compactEdgeScores();

    Graph G2(n);

    G2.indexEdges();

    G2.addEdge(0, 2);
    G2.addEdge(1, 2);
    G2.addEdge(2, 3);
    G2.addEdge(3, 5);
    G2.addEdge(4, 5);

    Betweenness centrality2(G2, false, true);
    centrality2.run();
    std::vector<double> bc2 = centrality2.edgeScores();
    for (size_t i = 0; i < bc.size(); i++) {
        EXPECT_EQ(bc[i], bc2[i]);
    }
}

TEST_F(BetweennessGTest, runApproxBetweennessSmallGraph) {
    /* Graph:
     0    3
      \  / \
       2   5
      / \ /
     1   4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    double epsilon = 0.1; // error
    double delta = 0.1;   // confidence
    ApproxBetweenness centrality(G, epsilon, delta);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    ASSERT_LE(centrality.scores().size(), 1.0 / (epsilon * epsilon));

    DEBUG("scores: ", bc);
}

TEST_F(BetweennessGTest, runApproxBetweenness) {
    DorogovtsevMendesGenerator generator(100);
    Graph G1 = generator.generate();
    Graph G(G1, true, false);
    ApproxBetweenness bc(G, 0.1, 0.1);
    bc.run();
    ApproxBetweenness bc1(G1, 0.1, 0.1);
    bc1.run();
}

TEST_F(BetweennessGTest, testBetweennessCentralityWeighted) {
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

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    const double tol = 1e-3;
    EXPECT_NEAR(0.0, bc[0], tol);
    EXPECT_NEAR(0.0, bc[1], tol);
    EXPECT_NEAR(23.0, bc[2], tol);
    EXPECT_NEAR(0.0, bc[3], tol);
    EXPECT_NEAR(0.0, bc[4], tol);
    EXPECT_NEAR(23.0, bc[5], tol);
    EXPECT_NEAR(0.0, bc[6], tol);
    EXPECT_NEAR(0.0, bc[7], tol);
}

TEST_F(BetweennessGTest, benchSequentialBetweennessCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    Betweenness bc(G);
    bc.run();
    std::vector<std::pair<node, double>> ranking = bc.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(BetweennessGTest, benchParallelBetweennessCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    Betweenness bc(G);
    bc.run();
    std::vector<std::pair<node, double>> ranking = bc.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(BetweennessGTest, runEstimateBetweenness) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");

    EstimateBetweenness abc2(G, 100);
    abc2.run();

    DEBUG("approximated betweenness scores: ", abc2.scores());
}

TEST_F(BetweennessGTest, testBetweennessMaximum) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    Graph G(6);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.indexEdges();

    Betweenness bet(G, false, false);
    bet.run();
    EXPECT_EQ(bet.maximum(), 10);
    EXPECT_NEAR(bet.centralization(), 1.79, 0.01);
}

TEST_F(BetweennessGTest, testBetweennessMaximumNormalized) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    Graph G(6);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.indexEdges();

    Betweenness betNormalized(G, true, true);
    betNormalized.run();
    EXPECT_EQ(betNormalized.maximum(), 1);
    EXPECT_NEAR(betNormalized.centralization(), 0.69, 0.01);
}

TEST_F(BetweennessGTest, testEdgeBetweennessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.indexEdges();

    Betweenness centrality(G, false, true);
    centrality.run();
    std::vector<double> bc = centrality.edgeScores();

    const double tol = 1e-3;
    EXPECT_NEAR(10.0, bc[0], tol);
    EXPECT_NEAR(10.0, bc[1], tol);
    EXPECT_NEAR(10.0, bc[2], tol);
    EXPECT_NEAR(10.0, bc[3], tol);
    EXPECT_NEAR(6.0, bc[4], tol);
    EXPECT_NEAR(6.0, bc[5], tol);
}

TEST_F(BetweennessGTest, debugEdgeBetweennessCentrality) {
    auto path = "input/PGPgiantcompo.graph";
    METISGraphReader reader;
    Graph G = reader.read(path);
    G.indexEdges();

    Betweenness centrality(G, false, true);
    centrality.run();
    std::vector<double> bc = centrality.edgeScores();
}

} /* namespace NetworKit */

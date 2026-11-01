// networkit-format
/*
 * KadabraBetweennessGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>

#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/KadabraBetweenness.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class KadabraBetweennessGTest : public testing::Test {};

TEST_F(KadabraBetweennessGTest, testKadabraEmptyGraph) {
    Graph G(0);
    EXPECT_THROW(KadabraBetweenness(G, 0.1, 0.1), std::invalid_argument);
}

TEST_F(KadabraBetweennessGTest, testKadabraAllNodesRemoved) {
    Graph G(3);
    G.removeNode(0);
    G.removeNode(1);
    G.removeNode(2);
    ASSERT_TRUE(G.isEmpty());
    ASSERT_GT(G.upperNodeIdBound(), 0);
    EXPECT_THROW(KadabraBetweenness(G, 0.1, 0.1), std::invalid_argument);
}

/**
 * This test succeeds with the fixed random seed (42).
 * However, the Kadabra algorithm computes a correct epsilon-approximation of
 * the betweenness centrality score of all the nodes of the graph with high
 * probability. Thus, it is possible that, for a different random seed, this
 * test fails.
 */
TEST_F(KadabraBetweennessGTest, testKadabraAbsolute) {
    Aux::Random::setSeed(42, true);
    const count n = 10;
    Graph g = ErdosRenyiGenerator(n, 0.1).generate();

    const double delta = 0.1;
    const double epsilon = 0.01;
    KadabraBetweenness kadabra(g, epsilon, delta);
    kadabra.run();
    auto scores = kadabra.topkScoresList();
    auto nodes = kadabra.topkNodesList();

    Betweenness betweenness(g, true);
    betweenness.run();
    count maxErrors = (count)std::ceil(delta * (double)n);

    count errors = 0;
    for (count i = 0; i < n; ++i) {
        if (std::abs(scores[i] - betweenness.score(nodes[i])) > delta) {
            ++errors;
        }
    }

    EXPECT_TRUE(errors <= maxErrors);
}

/**
 * This test succeeds with the fixed random seed (42).
 * However, the Kadabra algorithm finds the top-k nodes with
 * highest betweenness centrality with high probability. Thus, it is possible
 * that, for a different random seed, this test fails.
 */

TEST_F(KadabraBetweennessGTest, testKadabraTopK) {
    Aux::Random::setSeed(42, true);
    const count n = 10;
    Graph g = ErdosRenyiGenerator(n, 0.1).generate();

    const double delta = 0.1;
    const double epsilon = 0.01;
    const count k = 3;
    KadabraBetweenness kadabra(g, epsilon, delta, k);
    kadabra.run();
    auto kadabraRanking = kadabra.ranking();

    Betweenness betweenness(g, true);
    betweenness.run();
    auto betwRanking = betweenness.ranking();
    bool correctRanking = true;
    for (count i = 0; i < k; ++i) {
        if (betwRanking[i].first != kadabraRanking[i].first) {
            correctRanking = false;
            int j = static_cast<int>(i) - 1;
            while (j >= 0 && betwRanking[j].second == betwRanking[i].second) {
                --j;
            }
            ++j;
            while (j < static_cast<int>(n) && betwRanking[j].second == betwRanking[i].second) {
                if (betwRanking[j].first == kadabraRanking[i].first) {
                    correctRanking = true;
                    break;
                }
                ++j;
            }
        }
    }
    EXPECT_TRUE(correctRanking);
}

TEST_F(KadabraBetweennessGTest, testKadabraAbsoluteDeterministic) {
    const index seed = 42;
    Aux::Random::setSeed(seed, true);
    const count n = 300;
    Graph g = ErdosRenyiGenerator(n, 0.1).generate();

    const double delta = 0.1;
    const double epsilon = 0.01;
    Aux::Random::setSeed(seed, true);
    KadabraBetweenness kadabra(g, epsilon, delta, true);
    kadabra.run();
    auto scores = kadabra.topkScoresList();
    auto nodes = kadabra.topkNodesList();

    Aux::Random::setSeed(seed, true);
    KadabraBetweenness kadabra2(g, epsilon, delta, true);
    kadabra2.run();
    auto scores2 = kadabra2.topkScoresList();
    auto nodes2 = kadabra2.topkNodesList();

    for (count i = 0; i < n; ++i) {
        EXPECT_EQ(scores[i], scores2[i]);
    }
}

} /* namespace NetworKit */

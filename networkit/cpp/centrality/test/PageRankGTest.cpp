// networkit-format
/*
 * PageRankGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <cmath>
#include <utility>
#include <vector>

#include <gmock/gmock-more-matchers.h>
#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/PageRank.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>

#include "CentralityTestFixture.hpp"

namespace NetworKit {

class PageRankGTest : public DirectedWeightedGTest {};

INSTANTIATE_TEST_SUITE_P(InstantiationName, PageRankGTest, directedWeightedParams());

TEST_P(PageRankGTest, testPageRank) {
    SNAPGraphReader reader(isDirected());
    auto G = reader.read("input/wiki-Vote.txt");

    auto doTest = [&G](PageRank::Norm norm) {
        PageRank pr(G);
        pr.norm = norm;
        pr.run();

        auto pr_ranking = pr.ranking();
        constexpr double eps = 1e-3;
        if (G.isDirected()) {
            EXPECT_EQ(pr_ranking[0].first, 326);
            EXPECT_NEAR(pr_ranking[0].second, 0.00460, eps);
        } else {
            EXPECT_EQ(pr_ranking[0].first, 699);
            EXPECT_NEAR(pr_ranking[0].second, 0.00432, eps);
        }

        constexpr count maxIterations = 2;
        pr.maxIterations = maxIterations;
        pr.run();
        EXPECT_LE(pr.numberOfIterations(), maxIterations);
    };

    doTest(PageRank::Norm::L1_NORM);
    doTest(PageRank::Norm::L2_NORM);
}

TEST_F(PageRankGTest, testPageRankOnEmptyGraph) {
    const Graph G(0);

    PageRank pr(G);
    pr.run();

    EXPECT_THAT(pr.scores(), testing::IsEmpty());
    EXPECT_THAT(pr.ranking(), testing::IsEmpty());
    EXPECT_EQ(pr.numberOfIterations(), 0);
    EXPECT_DOUBLE_EQ(pr.maximum(), 0.0);
}

TEST_P(PageRankGTest, testNormalizedPageRank) {
    /* Graph:
     0 <---> 1
     \-> 2 <-/
     3       4

     Node 0,1 have directed edges to each other. Both have also an
     edge to node 2. Node 3 and 4 are isolated from the rest (sinks).
     This example is taken from "Comparing Apples and Oranges:
     Normalized PageRank for Evolving Graphs" by Berberich et al.
    */
    count n = 5;
    Graph G(n, isWeighted(), isDirected());
    G.addEdge(0, 1);
    G.addEdge(1, 0);
    G.addEdge(0, 2);
    G.addEdge(1, 2);

    auto doTest = [&G](PageRank::Norm norm) {
        PageRank pr(G, 0.85, 1e-8, true, PageRank::SinkHandling::DISTRIBUTE_SINKS);
        pr.norm = norm;
        pr.run();

        auto pr_scores = pr.scores();
        const double tol = 2e-4;

        // Values should be the same as in the original paper
        if (G.isDirected()) {
            EXPECT_NEAR(pr_scores[0], 1.7391, tol);
            EXPECT_NEAR(pr_scores[1], 1.7391, tol);
            EXPECT_NEAR(pr_scores[2], 2.4781, tol);
            EXPECT_NEAR(pr_scores[3], 1.0, tol);
            EXPECT_NEAR(pr_scores[4], 1.0, tol);
        } else {
            EXPECT_NEAR(pr_scores[0], 7.4026, tol);
            EXPECT_NEAR(pr_scores[1], 7.4026, tol);
            EXPECT_NEAR(pr_scores[2], 5.1948, tol);
            EXPECT_NEAR(pr_scores[3], 1.0, tol);
            EXPECT_NEAR(pr_scores[4], 1.0, tol);
        }
    };

    doTest(PageRank::Norm::L1_NORM);
    doTest(PageRank::Norm::L2_NORM);
}

TEST_F(PageRankGTest, testPageRankCentrality) {
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

    double damp = 0.85;
    PageRank centrality(G, damp);
    centrality.run();
    std::vector<double> cen = centrality.scores();

    // compare to Matlab results
    const double tol = 1e-4;
    EXPECT_NEAR(0.0753, std::fabs(cen[0]), tol);
    EXPECT_NEAR(0.0565, std::fabs(cen[1]), tol);
    EXPECT_NEAR(0.2552, std::fabs(cen[2]), tol);
    EXPECT_NEAR(0.1319, std::fabs(cen[3]), tol);
    EXPECT_NEAR(0.0942, std::fabs(cen[4]), tol);
    EXPECT_NEAR(0.2552, std::fabs(cen[5]), tol);
    EXPECT_NEAR(0.0753, std::fabs(cen[6]), tol);
    EXPECT_NEAR(0.0565, std::fabs(cen[7]), tol);
}

TEST_F(PageRankGTest, benchPageRankCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    double damp = 0.85;
    PageRank cen(G, damp);
    cen.run();
    std::vector<std::pair<node, double>> ranking = cen.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(PageRankGTest, benchNormalizedPageRankCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    double damp = 0.85;
    PageRank cen(G, damp, 1e-8, true);
    cen.run();
    std::vector<std::pair<node, double>> ranking = cen.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

} /* namespace NetworKit */

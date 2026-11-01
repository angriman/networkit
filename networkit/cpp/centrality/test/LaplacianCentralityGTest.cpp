// networkit-format
/*
 * LaplacianCentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <vector>

#include <gtest/gtest.h>

#include <networkit/centrality/LaplacianCentrality.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class LaplacianCentralityGTest : public testing::Test {};

TEST_F(LaplacianCentralityGTest, testLaplacianCentrality) {
    // The graph structure and reference values for the scores are taken from
    // Qi et al., Laplacian centrality: A new centrality measure for weighted
    // networks.
    //
    // See
    // https://math.wvu.edu/~cqzhang/Publication-files/my-paper/INS-2012-Laplacian-W.pdf.
    Graph G(6, true);

    G.addEdge(0, 1, 4);
    G.addEdge(0, 2, 2);
    G.addEdge(1, 2, 1);
    G.addEdge(1, 3, 2);
    G.addEdge(1, 4, 2);
    G.addEdge(4, 5, 1);

    LaplacianCentrality lc(G);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_EQ(140, scores[0]);
    EXPECT_EQ(180, scores[1]);
    EXPECT_EQ(56, scores[2]);
    EXPECT_EQ(44, scores[3]);
    EXPECT_EQ(52, scores[4]);
    EXPECT_EQ(8, scores[5]);
}

TEST_F(LaplacianCentralityGTest, testLaplacianCentralityNormalized) {
    Graph G(6, true);

    G.addEdge(0, 1, 4);
    G.addEdge(0, 2, 2);
    G.addEdge(1, 2, 1);
    G.addEdge(1, 3, 2);
    G.addEdge(1, 4, 2);
    G.addEdge(4, 5, 1);

    LaplacianCentrality lc(G, true);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_EQ(0.70, scores[0]);
    EXPECT_EQ(0.90, scores[1]);
    EXPECT_EQ(0.28, scores[2]);
    EXPECT_EQ(0.22, scores[3]);
    EXPECT_EQ(0.26, scores[4]);
    EXPECT_EQ(0.04, scores[5]);
}

TEST_F(LaplacianCentralityGTest, testLaplacianCentralityUnweighted) {
    Graph G(6);

    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(1, 4);
    G.addEdge(4, 5);

    LaplacianCentrality lc(G);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_EQ(18, scores[0]);
    EXPECT_EQ(34, scores[1]);
    EXPECT_EQ(18, scores[2]);
    EXPECT_EQ(10, scores[3]);
    EXPECT_EQ(16, scores[4]);
    EXPECT_EQ(6, scores[5]);
}

TEST_F(LaplacianCentralityGTest, testLaplacianCentralityFractionalWeights) {
    // The node's own weighted degree must not be truncated before being squared.
    // The weighted degrees here are 2.5, 3.75 and 1.25, and every value below is
    // exact in binary floating point.
    Graph G(3, true);

    G.addEdge(0, 1, 2.5);
    G.addEdge(1, 2, 1.25);

    LaplacianCentrality lc(G);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_DOUBLE_EQ(31.25, scores[0]);
    EXPECT_DOUBLE_EQ(37.5, scores[1]);
    EXPECT_DOUBLE_EQ(12.5, scores[2]);
}

TEST_F(LaplacianCentralityGTest, testLaplacianCentralityFractionalWeightsNormalized) {
    Graph G(3, true);

    G.addEdge(0, 1, 2.5);
    G.addEdge(1, 2, 1.25);

    LaplacianCentrality lc(G, true);
    lc.run();
    std::vector<double> scores = lc.scores();

    // Removing the middle vertex releases the whole Laplacian energy of this
    // path, so its normalized score is 1. Truncating the weighted degree shrank
    // the denominator and pushed the score above 1.
    EXPECT_NEAR(31.25 / 37.5, scores[0], 1e-12);
    EXPECT_NEAR(1.0, scores[1], 1e-12);
    EXPECT_NEAR(12.5 / 37.5, scores[2], 1e-12);
}

} /* namespace NetworKit */

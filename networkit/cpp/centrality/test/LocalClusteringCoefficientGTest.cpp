// networkit-format
/*
 * LocalClusteringCoefficientGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <vector>

#include <gtest/gtest.h>

#include <networkit/centrality/LocalClusteringCoefficient.hpp>
#include <networkit/centrality/LocalSquareClusteringCoefficient.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class LocalClusteringCoefficientGTest : public testing::Test {};

TEST_F(LocalClusteringCoefficientGTest, testLocalClusteringCoefficientUndirected) {
    count n = 16;
    Graph G(n, false, false);

    // 	// create graph used in Baur et al. and network analysis lecture
    G.addEdge(2, 4);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 7);
    G.addEdge(6, 7);

    G.addEdge(6, 8);
    G.addEdge(6, 9);
    G.addEdge(6, 11);
    G.addEdge(7, 12);
    G.addEdge(8, 9);

    G.addEdge(8, 10);
    G.addEdge(8, 11);
    G.addEdge(8, 13);
    G.addEdge(9, 10);
    G.addEdge(9, 11);

    G.addEdge(9, 13);
    G.addEdge(10, 11);
    G.addEdge(10, 13);
    G.addEdge(10, 14);
    G.addEdge(11, 13);

    G.addEdge(11, 14);
    G.addEdge(12, 15);
    G.addEdge(13, 14);
    G.addEdge(14, 15);

    EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
    EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

    // compute core decomposition
    LocalClusteringCoefficient lcc(G);
    lcc.run();
    std::vector<double> lccScores = lcc.scores();
    std::vector<double> reference = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                     0.5, 0.0, 0.8, 0.8, 0.8, 0.6666666666666666,
                                     0.0, 0.8, 0.5, 0.0};

    EXPECT_EQ(reference, lccScores);

    LocalClusteringCoefficient lccTurbo(G, true);
    lccTurbo.run();
    EXPECT_EQ(reference, lccTurbo.scores());

    // test throw runtime error for self-loop in graph
    Graph H(2);
    H.addEdge(0, 1);
    H.addEdge(1, 1);
    EXPECT_ANY_THROW(LocalClusteringCoefficient lcc(H));
}

TEST_F(LocalClusteringCoefficientGTest, testLocalClusteringCoefficientUndirected2) {
    Graph G(6, false, false);
    G.addEdge(1, 0);
    G.addEdge(2, 0);
    G.addEdge(2, 1);
    G.addEdge(3, 2);
    G.addEdge(3, 0);
    G.addEdge(3, 1);
    G.addEdge(4, 2);
    G.addEdge(4, 0);
    G.addEdge(5, 3);
    G.addEdge(5, 4);
    G.addEdge(5, 1);
    LocalClusteringCoefficient lcc(G);
    lcc.run();
    std::vector<double> lccScores = lcc.scores();
    std::vector<double> reference = {0.6666666666666666, 0.6666666666666666, 0.6666666666666666,
                                     0.6666666666666666, 0.3333333333333333, 0.3333333333333333};

    EXPECT_EQ(reference, lccScores);
}

TEST_F(LocalClusteringCoefficientGTest, testLocalSquareClusteringCoefficientUndirected) {
    Graph G(7, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 0);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(6, 3);
    LocalSquareClusteringCoefficient lscc(G);
    lscc.run();
    std::vector<double> lccScores = lscc.scores();
    std::vector<double> reference = {0.3333333333333333, 1.0, 0.3333333333333333, 0.2,
                                     0.3333333333333333, 1.0, 0.3333333333333333};

    EXPECT_EQ(reference, lccScores);
}

} /* namespace NetworKit */

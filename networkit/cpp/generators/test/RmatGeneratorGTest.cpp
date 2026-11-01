// networkit-format
/*
 * RmatGeneratorGTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <stdexcept>

#include <gtest/gtest.h>

#include <networkit/community/Modularity.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/generators/RmatGenerator.hpp>
#include <networkit/global/ClusteringCoefficient.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class RmatGeneratorGTest : public testing::Test {};

TEST_F(RmatGeneratorGTest, testRmatGeneratorException) {
    count scale = 9;
    count edgeFactor = 12;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.2;

    EXPECT_THROW(RmatGenerator rmat(scale, edgeFactor, a, b, c, d), std::runtime_error);
}

TEST_F(RmatGeneratorGTest, testRmatGenerator) {
    count scale = 9;
    count n = (1 << scale);
    count edgeFactor = 12;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.25;

    RmatGenerator rmat(scale, edgeFactor, a, b, c, d);
    Graph G = rmat.generate();

    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_LE(G.numberOfEdges(), n * edgeFactor);

    ClusteringCoefficient cc;
    double ccex = cc.exactGlobal(G);
    EXPECT_LE(ccex, 0.4);

    PLM clusterer(G, true);
    clusterer.run();
    Partition zeta = clusterer.getPartition();
    Modularity mod;
    double modVal = mod.getQuality(zeta, G);
    INFO("Modularity of R-MAT graph clustering: ", modVal);
    EXPECT_GE(modVal, 0.0);
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(RmatGeneratorGTest, testRmatGeneratorDistribution) {
    constexpr count scale = 2;
    constexpr count n = (1 << scale);
    count edgeFactor = 1;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.25;
    double edgeExpectedProbability[n][n] = {
        {0.0e+00, 0.0e+00, 0.0e+00, 0.0e+00},
        {2.242834007823646e-01, 0.0e+00, 0.0e+00, 0.0e+00},
        {2.242834007823646e-01, 1.0219127659785153e-01, 0.0e+00, 0.0e+00},
        {1.0219127659785153e-01, 1.7352532261978387e-01, 1.7352532261978387e-01, 0.0e+00},
    };

    RmatGenerator rmat(scale, edgeFactor, a, b, c, d, false, 0, true);
    count edgeCount[n][n]{{0}};
    count totalEdges = 0;
    // Now we generate a bunch of graphs and count the edges.
    for (index k = 0; k < 1000; k++) {
        Graph G = rmat.generate();
        G.forEdges([&edgeCount, &totalEdges](node u, node v) {
            edgeCount[u][v] += 1;
            totalEdges += 1;
        });
    }
    for (index i = 0; i < n; ++i) {
        for (index j = 0; j < n; ++j) {
            EXPECT_NEAR((static_cast<double>(edgeCount[i][j]) / static_cast<double>(totalEdges)),
                        edgeExpectedProbability[i][j], 0.01);
        }
    }
}

TEST_F(RmatGeneratorGTest, testRmatGeneratorReduceNodes) {
    count scale = 9;
    count n = (1 << scale);
    count edgeFactor = 12;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.25;
    int reducedNodes = 4;

    RmatGenerator rmat(scale, edgeFactor, a, b, c, d, false, reducedNodes);
    Graph G = rmat.generate();

    EXPECT_EQ(G.numberOfNodes(), n - reducedNodes);
    EXPECT_LE(G.numberOfEdges(), n * edgeFactor);

    ClusteringCoefficient cc;
    double ccex = cc.exactGlobal(G);
    EXPECT_LE(ccex, 0.4);

    PLM clusterer(G, true);
    clusterer.run();
    Partition zeta = clusterer.getPartition();
    Modularity mod;
    double modVal = mod.getQuality(zeta, G);
    INFO("Modularity of R-MAT graph clustering: ", modVal);
    EXPECT_GE(modVal, 0.0);
    EXPECT_TRUE(G.checkConsistency());
}

} /* namespace NetworKit */

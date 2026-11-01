// networkit-format
#include <cmath>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/centrality/HITS.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class HITSGTest : public testing::Test {};

TEST_F(HITSGTest, testDirectedToyGraph) {
    // Authorities are the principal eigenvector of [[2, 1], [1, 1]], hubs the one of
    // [[1, 1], [1, 2]]: both are (phi, 1) normalized, with phi being the golden ratio.
    Graph G(4, false, true);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(1, 3);

    HITS hits(G);
    hits.run();

    const double large = (1. + std::sqrt(5.)) / std::sqrt(10. + 2. * std::sqrt(5.));
    const double small = 2. / std::sqrt(10. + 2. * std::sqrt(5.));

    const auto &authorities = hits.authorityScores();
    EXPECT_NEAR(authorities[0], 0., 1e-6);
    EXPECT_NEAR(authorities[1], 0., 1e-6);
    EXPECT_NEAR(authorities[2], large, 1e-6);
    EXPECT_NEAR(authorities[3], small, 1e-6);

    const auto &hubs = hits.hubScores();
    EXPECT_NEAR(hubs[0], small, 1e-6);
    EXPECT_NEAR(hubs[1], large, 1e-6);
    EXPECT_NEAR(hubs[2], 0., 1e-6);
    EXPECT_NEAR(hubs[3], 0., 1e-6);

    EXPECT_EQ(hits.scores(), hits.authorityScores());
    EXPECT_GT(hits.numberOfIterations(), 1);
}

TEST_F(HITSGTest, testWeightsScaleHubScores) {
    Graph G(3, true, true);
    G.addEdge(0, 1, 2.);
    G.addEdge(2, 1, 1.);

    HITS hits(G);
    hits.run();

    EXPECT_NEAR(hits.scores()[1], 1., 1e-6);
    EXPECT_NEAR(hits.hubScores()[0], 2. / std::sqrt(5.), 1e-6);
    EXPECT_NEAR(hits.hubScores()[2], 1. / std::sqrt(5.), 1e-6);
}

TEST_F(HITSGTest, testUndirectedNonBipartiteScoresCoincide) {
    // Triangle with a pendant node.
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 0);
    G.addEdge(2, 3);

    HITS hits(G);
    hits.run();

    const auto &scores = hits.scores();
    G.forNodes([&](node u) -> void { EXPECT_NEAR(scores[u], hits.hubScores()[u], 1e-6); });

    // The scores are an eigenvector of the adjacency matrix.
    const double eigenvalue = scores[0] + scores[1] + scores[3];
    G.forNodes([&](node u) -> void {
        double neighborSum = 0.;
        G.forNeighborsOf(u, [&](node v) -> void { neighborSum += scores[v]; });
        EXPECT_NEAR(neighborSum, eigenvalue * scores[u] / scores[2], 1e-5);
    });
}

TEST_F(HITSGTest, testUndirectedBipartiteScoresDiffer) {
    // The largest eigenvalue of the squared adjacency matrix of a path is not simple, hence the
    // limit depends on the initial uniform scores and hubs differ from authorities.
    Graph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    HITS hits(G);
    hits.run();

    const std::vector<double> authorities{1. / std::sqrt(6.), 2. / std::sqrt(6.),
                                          1. / std::sqrt(6.)};
    for (node u = 0; u < 3; ++u) {
        EXPECT_NEAR(hits.scores()[u], authorities[u], 1e-6);
        EXPECT_NEAR(hits.hubScores()[u], 1. / std::sqrt(3.), 1e-6);
    }
}

TEST_F(HITSGTest, testScoresHaveUnitNorm) {
    Graph G(5, false, true);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 0);
    G.addEdge(4, 2);

    HITS hits(G);
    hits.run();

    double authorityNorm = 0., hubNorm = 0.;
    G.forNodes([&](node u) -> void {
        authorityNorm += hits.scores()[u] * hits.scores()[u];
        hubNorm += hits.hubScores()[u] * hits.hubScores()[u];
    });
    EXPECT_NEAR(authorityNorm, 1., 1e-9);
    EXPECT_NEAR(hubNorm, 1., 1e-9);
}

TEST_F(HITSGTest, testDeletedNodes) {
    Graph G(4, false, true);
    G.addEdge(0, 3);
    G.addEdge(1, 3);
    G.removeNode(1);

    HITS hits(G);
    hits.run();

    ASSERT_EQ(hits.scores().size(), 4);
    EXPECT_NEAR(hits.scores()[3], 1., 1e-6);
    EXPECT_NEAR(hits.hubScores()[0], 1., 1e-6);
    EXPECT_DOUBLE_EQ(hits.scores()[1], 0.);
    EXPECT_DOUBLE_EQ(hits.hubScores()[1], 0.);
}

TEST_F(HITSGTest, testEmptyGraph) {
    Graph G;
    HITS hits(G);
    hits.run();
    EXPECT_TRUE(hits.scores().empty());
    EXPECT_TRUE(hits.hubScores().empty());
    EXPECT_EQ(hits.numberOfIterations(), 0);
}

TEST_F(HITSGTest, testGraphWithoutEdges) {
    Graph G(3, false, true);
    HITS hits(G);
    hits.run();
    for (node u = 0; u < 3; ++u) {
        EXPECT_DOUBLE_EQ(hits.scores()[u], 0.);
        EXPECT_DOUBLE_EQ(hits.hubScores()[u], 0.);
    }
}

TEST_F(HITSGTest, testMaxIterations) {
    Graph G(4, false, true);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(1, 3);

    HITS hits(G, 1e-12, 2);
    hits.run();
    EXPECT_EQ(hits.numberOfIterations(), 2);
}

TEST_F(HITSGTest, testRunTwiceGivesSameScores) {
    Graph G(3, false, true);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(0, 2);

    HITS hits(G);
    hits.run();
    const auto authorities = hits.scores();
    const auto hubs = hits.hubScores();
    hits.run();
    EXPECT_EQ(hits.scores(), authorities);
    EXPECT_EQ(hits.hubScores(), hubs);
}

TEST_F(HITSGTest, testInvalidParameters) {
    Graph G(2);
    EXPECT_THROW(HITS(G, 0.), std::invalid_argument);
    EXPECT_THROW(HITS(G, 1e-8, 0), std::invalid_argument);
}

TEST_F(HITSGTest, testScoresBeforeRunThrow) {
    Graph G(2);
    HITS hits(G);
    EXPECT_THROW(hits.hubScores(), std::runtime_error);
    EXPECT_THROW(hits.numberOfIterations(), std::runtime_error);
}

} /* namespace NetworKit */

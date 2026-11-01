// networkit-format
/*
 * ClosenessGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <limits>
#include <vector>

#include <gmock/gmock-more-matchers.h>
#include <gtest/gtest.h>

#include <networkit/centrality/ApproxCloseness.hpp>
#include <networkit/centrality/ApproxElectricalCloseness.hpp>
#include <networkit/centrality/Closeness.hpp>
#include <networkit/centrality/DynTopHarmonicCloseness.hpp>
#include <networkit/centrality/HarmonicCloseness.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>

#include "CentralityTestFixture.hpp"

namespace NetworKit {

class ClosenessGTest : public DirectedWeightedGTest {};

INSTANTIATE_TEST_SUITE_P(InstantiationName, ClosenessGTest, directedWeightedParams());

// FIXME look out for tolerance limit in paper sample nodes
TEST_F(ClosenessGTest, testApproxClosenessCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");

    ApproxCloseness acc(G, 453, 0, true);
    acc.run();

    std::vector<double> acc_scores = acc.scores();

    ASSERT_EQ(acc_scores.size(), 453);

    // compare sampled approx closeness values vs real closeness values
    // its a good compromise to not let the exact closness algorithm run
    ASSERT_NEAR(acc_scores[0], 0.416206, 0.000001);
    ASSERT_NEAR(acc_scores[10], 0.355906, 0.000001);
    ASSERT_NEAR(acc_scores[87], 0.420465, 0.000001);
    ASSERT_NEAR(acc_scores[121], 0.38865, 0.000001);
    ASSERT_NEAR(acc_scores[178], 0.4, 0.000001);
    ASSERT_NEAR(acc_scores[254], 0.397188, 0.000001);
    ASSERT_NEAR(acc_scores[307], 0.398238, 0.000001);
    ASSERT_NEAR(acc_scores[398], 0.37604, 0.000001);
    ASSERT_NEAR(acc_scores[406], 0.360734, 0.000001);
    ASSERT_NEAR(acc_scores[446], 0.396491, 0.000001);
}

TEST_F(ClosenessGTest, testApproxClosenessCentralityOnToyGraph) {
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

    ApproxCloseness acc(G, 6, 0.1, false);
    acc.run();
    std::vector<double> cc = acc.scores();

    double maximum = acc.maximum();

    const double tol = 0.2;
    EXPECT_NEAR(0.1, cc[0], tol);
    EXPECT_NEAR(0.1, cc[1], tol);
    EXPECT_NEAR(0.166667, cc[2], tol);
    EXPECT_NEAR(0.125, cc[3], tol);
    EXPECT_NEAR(0.125, cc[4], tol);
    EXPECT_NEAR(0.1, cc[5], tol);
    EXPECT_NEAR(0.2, maximum, tol);

    ApproxCloseness acc2(G, 4, 0.1, true);
    acc2.run();
    std::vector<double> cc2 = acc2.scores();

    double maximum2 = acc2.maximum();

    EXPECT_NEAR(0.5, cc2[0], tol);
    EXPECT_NEAR(0.5, cc2[1], tol);
    EXPECT_NEAR(0.833335, cc2[2], tol);
    EXPECT_NEAR(0.625, cc2[3], tol);
    EXPECT_NEAR(0.625, cc2[4], tol);
    EXPECT_NEAR(0.5, cc2[5], tol);
    EXPECT_NEAR(0.2, maximum2, tol);
}

TEST_F(ClosenessGTest, testApproxClosenessCentralityOnDirectedToyGraph) {
    /* Graph:
     0 -> 3
     ^   /
      \ /
       2    5
      /  \ /
     1    4
    */
    Graph G(6, false, true);

    G.addEdge(0, 3);
    G.addEdge(2, 0);
    G.addEdge(1, 2);
    G.addEdge(2, 1);
    G.addEdge(3, 2);
    G.addEdge(2, 3);
    G.addEdge(4, 2);
    G.addEdge(2, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 4);

    ApproxCloseness acc(G, 6, 0.1, true);
    acc.run();
    std::vector<double> cc = acc.scores();
    double tol = 0.01;

    EXPECT_NEAR(acc.maximum(), 0.2, tol);
    EXPECT_NEAR(cc[0], 0.38, tol);
    EXPECT_NEAR(cc[1], 0.50, tol);
    EXPECT_NEAR(cc[2], 0.83, tol);
    EXPECT_NEAR(cc[3], 0.50, tol);
    EXPECT_NEAR(cc[4], 0.62, tol);
    EXPECT_NEAR(cc[5], 0.41, tol);
}

TEST_P(ClosenessGTest, testClosenessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n, isWeighted(), isDirected());

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    if (isDirected()) {
        G.addEdge(4, 1);
        G.addEdge(3, 0);
        G.addEdge(5, 2);
    }

    if (isWeighted()) {
        Aux::Random::setSeed(42, false);
        GraphTools::randomizeWeights(G);
    }

    auto computeCloseness = [&](node u, ClosenessVariant variant, bool normalized) {
        Dijkstra dij(G, u, false);
        dij.run();
        auto dists = dij.getDistances();

        double sumDist = 0.;
        double reached = 0.;

        G.forNodes([&](node v) {
            if (variant == ClosenessVariant::STANDARD) {
                sumDist += dists[v];
            } else if (dists[v] != std::numeric_limits<double>::max()) {
                ++reached;
                sumDist += dists[v];
            }
        });

        double score = sumDist;
        if (score) {
            score = (variant == ClosenessVariant::STANDARD)
                        ? 1.0 / score
                        : (reached - 1.0) / sumDist / (G.numberOfNodes() - 1.0);
        }

        if (normalized) {
            score *= (variant == ClosenessVariant::STANDARD ? G.numberOfNodes() : reached) - 1.0;
        }

        return score;
    };

    for (auto variant : {ClosenessVariant::STANDARD, ClosenessVariant::GENERALIZED}) {
        for (auto normalized : {true, false}) {
            Closeness centrality(G, normalized, variant);
            centrality.run();
            const auto bc = centrality.scores();

            G.forNodes(
                [&](node u) { EXPECT_DOUBLE_EQ(bc[u], computeCloseness(u, variant, normalized)); });
        }
    }
}

TEST_F(ClosenessGTest, testHarmonicClosenessSingleton) {
    for (bool normalized : {false, true}) {
        for (bool directed : {false, true}) {
            Graph G(1, false, directed);
            HarmonicCloseness centrality(G, normalized);
            centrality.run();
            EXPECT_DOUBLE_EQ(centrality.score(0), 0.);
            EXPECT_DOUBLE_EQ(centrality.maximum(), 0.);
        }
        Graph empty(0);
        HarmonicCloseness centrality(empty, normalized);
        centrality.run();
        EXPECT_THAT(centrality.scores(), testing::IsEmpty());
        EXPECT_DOUBLE_EQ(centrality.maximum(), 0.);
    }
}

TEST_F(ClosenessGTest, testHarmonicClosenessMaximum) {
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);

    for (bool normalized : {false, true}) {
        HarmonicCloseness centrality(G, normalized);
        centrality.run();
        EXPECT_DOUBLE_EQ(centrality.maximum(), normalized ? 1. : 3.);
        // Raw scores are [11/6, 5/2, 5/2, 11/6]. Scaling every score
        // and its theoretical maximum leaves centralization at 2/5.
        EXPECT_NEAR(centrality.centralization(), 2. / 5., 1e-12);
    }
}

TEST_F(ClosenessGTest, testHarmonicClosenessCentrality) {
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

    HarmonicCloseness centrality(G, false);
    centrality.run();
    std::vector<double> hc = centrality.scores();

    double maximum = centrality.maximum();

    const double tol = 1e-3;
    EXPECT_NEAR(2.833, hc[0], tol);
    EXPECT_NEAR(2.833, hc[1], tol);
    EXPECT_NEAR(4.5, hc[2], tol);
    EXPECT_NEAR(3.5, hc[3], tol);
    EXPECT_NEAR(3.5, hc[4], tol);
    EXPECT_NEAR(3.1667, hc[5], tol);
    EXPECT_NEAR(5., maximum, tol);
}

TEST_P(ClosenessGTest, testDynTopHarmonicCloseness) {
    auto G1 = DorogovtsevMendesGenerator(500).generate();
    Graph G(G1, false, isDirected());

    constexpr count k = 10;

    DynTopHarmonicCloseness centrality(G, k, false);
    centrality.run();

    HarmonicCloseness reference(G, false);
    reference.run();

    auto scores = centrality.ranking();
    auto refScores = reference.ranking();

    for (count j = 0; j < k; ++j) {
        EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }

    count numInsertions = 1;

    std::vector<GraphEvent> deletions;
    std::vector<GraphEvent> insertions;

    for (count i = 0; i < numInsertions; i++) {

        node u = G.upperNodeIdBound();
        node v = G.upperNodeIdBound();

        do {
            u = GraphTools::randomNode(G);
            v = GraphTools::randomNode(G);
        } while (G.hasEdge(u, v));

        GraphEvent edgeAddition(GraphEvent::EDGE_ADDITION, u, v);
        insertions.insert(insertions.begin(), edgeAddition);

        GraphEvent edgeDeletion(GraphEvent::EDGE_REMOVAL, u, v);
        deletions.push_back(edgeDeletion);

        G.addEdge(u, v);
    }

    for (auto &e : insertions) {
        G.removeEdge(e.u, e.v);
    }

    for (auto &edgeAddition : insertions) {

        node u = edgeAddition.u;
        node v = edgeAddition.v;

        G.addEdge(u, v);

        centrality.update(edgeAddition);
        reference.run();

        scores = centrality.ranking();
        refScores = reference.ranking();

        for (count j = 0; j < k; ++j) {
            EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
        }
    }

    for (auto &edgeDeletion : deletions) {
        node u = edgeDeletion.u;
        node v = edgeDeletion.v;

        G.removeEdge(u, v);

        centrality.update(edgeDeletion);
        reference.run();

        scores = centrality.ranking();
        refScores = reference.ranking();

        for (count j = 0; j < k; ++j) {
            EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
        }
    }

    for (auto &edgeInsertion : insertions) {
        G.addEdge(edgeInsertion.u, edgeInsertion.v);
    }

    reference.run();
    centrality.updateBatch(insertions);

    scores = centrality.ranking();
    refScores = reference.ranking();

    for (count j = 0; j < k; ++j) {
        EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }
}

TEST_F(ClosenessGTest, testApproxElectricalCloseness) {
    const double eps = 0.1;
    const count n = 75;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        auto G = HyperbolicGenerator(n, 10, 3).generate();
        G = ConnectedComponents::extractLargestConnectedComponent(G, true);
        auto new_n = G.numberOfNodes();

        // Create a biconnected component with size 2.
        G.addNodes(2);
        G.addEdge(new_n - 1, new_n);
        G.addEdge(new_n, new_n + 1);

        ApproxElectricalCloseness apx(G);
        apx.run();
        const auto diag = apx.getDiagonal();
        const auto gt = apx.computeExactDiagonal(1e-12);
        G.forNodes([&](node u) { EXPECT_NEAR(diag[u], gt[u], eps); });
        EXPECT_EQ(apx.scores().size(), G.numberOfNodes());
    }
}

} /* namespace NetworKit */

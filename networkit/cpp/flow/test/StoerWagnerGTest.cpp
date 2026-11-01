// networkit-format
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/flow/StoerWagner.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class StoerWagnerGTest : public testing::Test {
protected:
    // Total weight of the edges with exactly one end point in side.
    static edgeweight cutWeight(const Graph &G, const std::vector<node> &side) {
        std::vector<bool> inSide(G.upperNodeIdBound(), false);
        for (const node u : side) {
            inSide[u] = true;
        }
        edgeweight weight = 0;
        G.forEdges([&](node u, node v, edgeweight w) -> void {
            if (inSide[u] != inSide[v]) {
                weight += w;
            }
        });
        return weight;
    }

    // Minimum over all bipartitions of the nodes of G, which must be consecutive.
    static edgeweight bruteForceMinCut(const Graph &G) {
        const count n = G.numberOfNodes();
        edgeweight best = std::numeric_limits<edgeweight>::infinity();
        // Node n - 1 is always on the other side, so each bipartition is enumerated once.
        for (count mask = 1; mask < (count{1} << (n - 1)); ++mask) {
            std::vector<node> side;
            for (node u = 0; u + 1 < n; ++u) {
                if (mask & (count{1} << u)) {
                    side.push_back(u);
                }
            }
            best = std::min(best, cutWeight(G, side));
        }
        return best;
    }
};

TEST_F(StoerWagnerGTest, testExampleFromPaper) {
    Graph G(8, true);
    G.addEdge(0, 1, 2);
    G.addEdge(0, 4, 3);
    G.addEdge(1, 2, 3);
    G.addEdge(1, 4, 2);
    G.addEdge(1, 5, 2);
    G.addEdge(2, 3, 4);
    G.addEdge(2, 6, 2);
    G.addEdge(3, 6, 2);
    G.addEdge(3, 7, 2);
    G.addEdge(4, 5, 3);
    G.addEdge(5, 6, 1);
    G.addEdge(6, 7, 3);

    StoerWagner algo(G);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getCutValue(), 4);
    const auto &side = algo.getCutSide();
    EXPECT_TRUE(side == std::vector<node>({2, 3, 6, 7}) || side == std::vector<node>({0, 1, 4, 5}));
}

TEST_F(StoerWagnerGTest, testUnweightedBridge) {
    // Two triangles connected by a single edge.
    Graph G(6);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 0);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 3);
    G.addEdge(2, 3);

    StoerWagner algo(G);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getCutValue(), 1);
    EXPECT_EQ(algo.getCutSide().size(), 3);
    EXPECT_DOUBLE_EQ(cutWeight(G, algo.getCutSide()), 1);
}

TEST_F(StoerWagnerGTest, testTwoNodes) {
    Graph G(2, true);
    G.addEdge(0, 1, 2.5);

    StoerWagner algo(G);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getCutValue(), 2.5);
    EXPECT_EQ(algo.getCutSide().size(), 1);
}

TEST_F(StoerWagnerGTest, testDisconnectedGraph) {
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(2, 3);

    StoerWagner algo(G);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getCutValue(), 0);
    EXPECT_DOUBLE_EQ(cutWeight(G, algo.getCutSide()), 0);
    EXPECT_FALSE(algo.getCutSide().empty());
    EXPECT_LT(algo.getCutSide().size(), 4);
}

TEST_F(StoerWagnerGTest, testSelfLoopsAndParallelEdges) {
    Graph G(3, true);
    G.addEdge(0, 0, 100);
    G.addEdge(0, 1, 1);
    G.addEdge(0, 1, 1);
    G.addEdge(1, 2, 3);

    StoerWagner algo(G);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getCutValue(), 2);
    EXPECT_EQ(algo.getCutSide(), std::vector<node>({0}));
}

TEST_F(StoerWagnerGTest, testDeletedNodes) {
    Graph G(5);
    G.addEdge(0, 2);
    G.addEdge(2, 4);
    G.addEdge(4, 0);
    G.addEdge(4, 3);
    G.addEdge(1, 2);
    G.removeNode(1);

    StoerWagner algo(G);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getCutValue(), 1);
    EXPECT_EQ(algo.getCutSide(), std::vector<node>({3}));
}

TEST_F(StoerWagnerGTest, testMatchesBruteForce) {
    for (const count seed : {1, 2, 3, 4, 5, 6, 7, 8}) {
        Aux::Random::setSeed(seed, false);
        for (const double probability : {0.3, 0.6, 0.9}) {
            Graph G = ErdosRenyiGenerator(10, probability).generate();
            G = GraphTools::toWeighted(G);
            G.forEdges([&](node u, node v) -> void {
                G.setWeight(u, v, static_cast<edgeweight>(Aux::Random::integer(1, 9)));
            });

            StoerWagner algo(G);
            algo.run();

            EXPECT_DOUBLE_EQ(algo.getCutValue(), bruteForceMinCut(G));
            EXPECT_DOUBLE_EQ(algo.getCutValue(), cutWeight(G, algo.getCutSide()));
            EXPECT_FALSE(algo.getCutSide().empty());
            EXPECT_LT(algo.getCutSide().size(), G.numberOfNodes());
        }
    }
}

TEST_F(StoerWagnerGTest, testRunTwice) {
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);

    StoerWagner algo(G);
    algo.run();
    const auto side = algo.getCutSide();
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getCutValue(), 1);
    EXPECT_EQ(algo.getCutSide(), side);
}

TEST_F(StoerWagnerGTest, testInvalidInput) {
    EXPECT_THROW(StoerWagner{Graph(3, false, true)}, std::runtime_error);
    EXPECT_THROW(StoerWagner{Graph(1)}, std::runtime_error);

    Graph negative(2, true);
    negative.addEdge(0, 1, -1);
    EXPECT_THROW(StoerWagner{negative}, std::runtime_error);
}

TEST_F(StoerWagnerGTest, testResultsBeforeRunThrow) {
    Graph G(2);
    StoerWagner algo(G);
    EXPECT_THROW(algo.getCutValue(), std::runtime_error);
    EXPECT_THROW(algo.getCutSide(), std::runtime_error);
}

} // namespace NetworKit

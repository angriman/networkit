// networkit-format
#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/SteinerTree.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace NetworKit {

class SteinerTreeGTest : public testing::Test {
protected:
    static void expectValidTree(const Graph &G, const Graph &tree,
                                const std::vector<node> &terminals) {
        UnionFind components(G.upperNodeIdBound());
        count nodesInTree = 0;
        tree.forNodes([&](node u) -> void {
            if (tree.degree(u) > 0) {
                ++nodesInTree;
            }
            // Non-terminal leaves would only add weight.
            if (tree.degree(u) == 1) {
                EXPECT_NE(std::find(terminals.begin(), terminals.end(), u), terminals.end());
            }
        });
        tree.forEdges([&](node u, node v, edgeweight weight) -> void {
            EXPECT_TRUE(G.hasEdge(u, v));
            EXPECT_DOUBLE_EQ(G.weight(u, v), weight);
            EXPECT_NE(components.find(u), components.find(v)) << "cycle through " << u;
            components.merge(u, v);
        });
        if (terminals.size() > 1) {
            EXPECT_EQ(tree.numberOfEdges() + 1, nodesInTree);
        }
        for (const node terminal : terminals) {
            EXPECT_EQ(components.find(terminal), components.find(terminals.front()));
        }
    }

    // Lightest spanning tree over all node subsets that contain the terminals. Node ids must be
    // consecutive and small.
    static edgeweight optimalWeight(const Graph &G, const std::vector<node> &terminals) {
        const count n = G.numberOfNodes();
        std::vector<std::tuple<edgeweight, node, node>> edges;
        G.forEdges(
            [&](node u, node v, edgeweight weight) -> void { edges.emplace_back(weight, u, v); });
        std::sort(edges.begin(), edges.end());

        uint64_t requiredNodes = 0;
        for (const node terminal : terminals) {
            requiredNodes |= uint64_t{1} << terminal;
        }

        edgeweight best = std::numeric_limits<edgeweight>::infinity();
        for (uint64_t subset = 0; subset < (uint64_t{1} << n); ++subset) {
            if ((subset & requiredNodes) != requiredNodes) {
                continue;
            }
            UnionFind components(n);
            edgeweight weight = 0;
            count remaining = 0;
            for (node u = 0; u < n; ++u) {
                remaining += (subset >> u) & 1;
            }
            for (const auto &[w, u, v] : edges) {
                if (((subset >> u) & 1) && ((subset >> v) & 1)
                    && components.find(u) != components.find(v)) {
                    components.merge(u, v);
                    weight += w;
                    --remaining;
                }
            }
            if (remaining == 1) {
                best = std::min(best, weight);
            }
        }
        return best;
    }
};

TEST_F(SteinerTreeGTest, testUsesASteinerNode) {
    // A star is lighter than any tree on the terminals alone.
    Graph G(4, true);
    G.addEdge(0, 3, 1);
    G.addEdge(1, 3, 1);
    G.addEdge(2, 3, 1);
    G.addEdge(0, 1, 2.5);
    G.addEdge(1, 2, 2.5);
    const std::vector<node> terminals{0, 1, 2};

    SteinerTree algo(G, terminals);
    algo.run();

    ASSERT_TRUE(algo.hasTree());
    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 3);
    EXPECT_EQ(algo.getTree().degree(3), 3);
    expectValidTree(G, algo.getTree(), terminals);
}

TEST_F(SteinerTreeGTest, testPathBetweenTwoTerminalsIsShortest) {
    Graph G(5, true);
    G.addEdge(0, 1, 1);
    G.addEdge(1, 2, 1);
    G.addEdge(2, 3, 1);
    G.addEdge(0, 4, 2);
    G.addEdge(4, 3, 5);
    const std::vector<node> terminals{0, 3};

    SteinerTree algo(G, terminals);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 3);
    expectValidTree(G, algo.getTree(), terminals);
}

TEST_F(SteinerTreeGTest, testUnweightedGraph) {
    Graph G(5);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(3, 4);

    SteinerTree algo(G, {0, 2, 4});
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 4);
    EXPECT_TRUE(algo.getTree().isWeighted());
}

TEST_F(SteinerTreeGTest, testAllNodesAreTerminals) {
    // The result is a spanning tree, and the guarantee makes it at most twice as heavy as optimal.
    Graph G(4, true);
    G.addEdge(0, 1, 1);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(3, 0, 4);
    const std::vector<node> terminals{0, 1, 2, 3};

    SteinerTree algo(G, terminals);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 6);
    expectValidTree(G, algo.getTree(), terminals);
}

TEST_F(SteinerTreeGTest, testDisconnectedTerminals) {
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(2, 3);

    SteinerTree algo(G, {0, 1, 2});
    algo.run();

    EXPECT_FALSE(algo.hasTree());
    EXPECT_EQ(algo.getTree().numberOfEdges(), 0);
    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 0);
}

TEST_F(SteinerTreeGTest, testComponentsWithoutTerminalsAreIgnored) {
    Graph G(5);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(3, 4);

    SteinerTree algo(G, {0, 2});
    algo.run();

    ASSERT_TRUE(algo.hasTree());
    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 2);
}

TEST_F(SteinerTreeGTest, testFewTerminals) {
    Graph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    for (const std::vector<node> &terminals :
         {std::vector<node>{}, std::vector<node>{1}, std::vector<node>{1, 1}}) {
        SteinerTree algo(G, terminals);
        algo.run();
        EXPECT_TRUE(algo.hasTree());
        EXPECT_EQ(algo.getTree().numberOfEdges(), 0);
        EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 0);
    }
}

TEST_F(SteinerTreeGTest, testZeroWeightsAndSelfLoops) {
    Graph G(4, true);
    G.addEdge(0, 1, 0);
    G.addEdge(1, 2, 0);
    G.addEdge(1, 1, 3);
    G.addEdge(2, 3, 1);
    const std::vector<node> terminals{0, 3};

    SteinerTree algo(G, terminals);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 1);
    expectValidTree(G, algo.getTree(), terminals);
}

TEST_F(SteinerTreeGTest, testDeletedNodes) {
    Graph G(5, true);
    G.addEdge(0, 2, 1);
    G.addEdge(2, 4, 1);
    G.addEdge(1, 0, 1);
    G.addEdge(1, 4, 1);
    G.removeNode(1);
    G.removeNode(3);

    SteinerTree algo(G, {0, 4});
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 2);
    EXPECT_FALSE(algo.getTree().hasNode(1));
}

TEST_F(SteinerTreeGTest, testApproximationGuarantee) {
    for (const count seed : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}) {
        Aux::Random::setSeed(seed, false);
        for (const count numberOfTerminals : {2, 3, 5}) {
            constexpr count n = 10;
            Graph G(n, true);
            for (node u = 1; u < n; ++u) {
                G.addEdge(Aux::Random::integer(u - 1), u,
                          static_cast<edgeweight>(Aux::Random::integer(1, 9)));
            }
            for (node u = 0; u < n; ++u) {
                for (node v = u + 1; v < n; ++v) {
                    if (!G.hasEdge(u, v) && Aux::Random::probability() < 0.3) {
                        G.addEdge(u, v, static_cast<edgeweight>(Aux::Random::integer(1, 9)));
                    }
                }
            }

            std::vector<node> nodes(n);
            std::iota(nodes.begin(), nodes.end(), node{0});
            std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());
            const std::vector<node> terminals(nodes.begin(), nodes.begin() + numberOfTerminals);

            SteinerTree algo(G, terminals);
            algo.run();

            ASSERT_TRUE(algo.hasTree());
            expectValidTree(G, algo.getTree(), terminals);
            const edgeweight optimal = optimalWeight(G, terminals);
            EXPECT_GE(algo.getTotalWeight(), optimal - 1e-9);
            const double bound = 2. * (1. - 1. / static_cast<double>(numberOfTerminals));
            EXPECT_LE(algo.getTotalWeight(), bound * optimal + 1e-9);
        }
    }
}

TEST_F(SteinerTreeGTest, testRunTwice) {
    Graph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    SteinerTree algo(G, {0, 2});
    algo.run();
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 2);
    EXPECT_EQ(algo.getTree().numberOfEdges(), 2);
}

TEST_F(SteinerTreeGTest, testInvalidInput) {
    EXPECT_THROW(SteinerTree(Graph(2, false, true), {0}), std::runtime_error);
    EXPECT_THROW(SteinerTree(Graph(2), {5}), std::runtime_error);

    Graph negative(2, true);
    negative.addEdge(0, 1, -1);
    EXPECT_THROW(SteinerTree(negative, {0}), std::runtime_error);
}

TEST_F(SteinerTreeGTest, testResultsBeforeRunThrow) {
    Graph G(2);
    SteinerTree algo(G, {0});
    EXPECT_THROW(algo.getTree(), std::runtime_error);
    EXPECT_THROW(algo.getTotalWeight(), std::runtime_error);
    EXPECT_THROW(algo.hasTree(), std::runtime_error);
}

} /* namespace NetworKit */

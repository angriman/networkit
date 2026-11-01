// networkit-format
#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/EulerianPath.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class EulerianPathGTest : public testing::TestWithParam<bool> {
protected:
    bool isDirected() const { return GetParam(); }

    // Checks that the path uses each edge exactly once.
    static void expectValidPath(const Graph &G, const std::vector<node> &path) {
        ASSERT_EQ(path.size(), G.numberOfEdges() + 1);

        auto key = [&](node u, node v) -> std::pair<node, node> {
            if (G.isDirected()) {
                return {u, v};
            }
            return {std::min(u, v), std::max(u, v)};
        };

        std::map<std::pair<node, node>, count> remaining;
        G.forEdges([&](node u, node v) -> void { ++remaining[key(u, v)]; });

        for (index i = 0; i + 1 < path.size(); ++i) {
            auto &multiplicity = remaining[key(path[i], path[i + 1])];
            ASSERT_GT(multiplicity, 0);
            --multiplicity;
        }
    }
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, EulerianPathGTest, testing::Bool());

TEST_P(EulerianPathGTest, testCircuit) {
    Graph G(4, false, isDirected());
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 0);

    EulerianPath algo(G);
    algo.run();

    ASSERT_TRUE(algo.hasEulerianPath());
    expectValidPath(G, algo.getPath());
    EXPECT_EQ(algo.getPath().front(), algo.getPath().back());
}

TEST_P(EulerianPathGTest, testOpenPathStartsAtTheRightNode) {
    // A cycle with a tail: the path has to start at the end of the tail.
    Graph G(5, false, isDirected());
    G.addEdge(4, 0);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 0);

    EulerianPath algo(G);
    algo.run();

    ASSERT_TRUE(algo.hasEulerianPath());
    expectValidPath(G, algo.getPath());
    if (isDirected()) {
        EXPECT_EQ(algo.getPath().front(), 4);
        EXPECT_EQ(algo.getPath().back(), 0);
    }
}

TEST_P(EulerianPathGTest, testNestedCyclesNeedSplicing) {
    // Two cycles through node 0: a greedy walk can close the first one too early.
    Graph G(5, false, isDirected());
    G.addEdge(0, 1);
    G.addEdge(1, 0);
    G.addEdge(0, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(4, 0);

    EulerianPath algo(G);
    algo.run();

    ASSERT_TRUE(algo.hasEulerianPath());
    expectValidPath(G, algo.getPath());
}

TEST_P(EulerianPathGTest, testSelfLoopsAndParallelEdges) {
    Graph G(3, false, isDirected());
    G.addEdge(0, 0);
    G.addEdge(0, 1);
    G.addEdge(0, 1);
    G.addEdge(1, 1);
    G.addEdge(1, 2);

    EulerianPath algo(G);
    algo.run();

    if (isDirected()) {
        EXPECT_FALSE(algo.hasEulerianPath());
    } else {
        ASSERT_TRUE(algo.hasEulerianPath());
        expectValidPath(G, algo.getPath());
    }
}

TEST_P(EulerianPathGTest, testDisconnectedEdgesHaveNoPath) {
    // The degrees of two disjoint cycles allow a circuit.
    Graph G(6, false, isDirected());
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 0);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 3);

    EulerianPath algo(G);
    algo.run();

    EXPECT_FALSE(algo.hasEulerianPath());
    EXPECT_TRUE(algo.getPath().empty());
}

TEST_P(EulerianPathGTest, testIsolatedAndDeletedNodesAreIgnored) {
    Graph G(6, false, isDirected());
    G.addEdge(1, 3);
    G.addEdge(3, 4);
    G.addEdge(4, 1);
    G.addEdge(0, 5);
    G.removeNode(0);

    EulerianPath algo(G);
    algo.run();

    ASSERT_TRUE(algo.hasEulerianPath());
    expectValidPath(G, algo.getPath());
}

TEST_P(EulerianPathGTest, testGraphWithoutEdges) {
    Graph G(3, false, isDirected());
    EulerianPath algo(G);
    algo.run();
    EXPECT_TRUE(algo.hasEulerianPath());
    EXPECT_TRUE(algo.getPath().empty());
}

TEST_P(EulerianPathGTest, testRandomWalkGraphs) {
    // The edges of a random walk trivially have an Eulerian path.
    for (const count seed : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}) {
        Aux::Random::setSeed(seed, false);
        constexpr count n = 30, steps = 300;
        Graph G(n, false, isDirected());
        node current = Aux::Random::index(n);
        for (count step = 0; step < steps; ++step) {
            const node next = Aux::Random::index(n);
            G.addEdge(current, next);
            current = next;
        }

        EulerianPath algo(G);
        algo.run();

        ASSERT_TRUE(algo.hasEulerianPath());
        expectValidPath(G, algo.getPath());
    }
}

TEST_F(EulerianPathGTest, testUndirectedOddNodes) {
    // Star with three leaves: four nodes of odd degree.
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(0, 3);

    EulerianPath algo(G);
    algo.run();

    EXPECT_FALSE(algo.hasEulerianPath());
}

TEST_F(EulerianPathGTest, testDirectedUnbalancedNodes) {
    // Node 0 has two more outgoing than incoming edges.
    Graph G(3, false, true);
    G.addEdge(0, 1);
    G.addEdge(0, 2);

    EulerianPath algo(G);
    algo.run();

    EXPECT_FALSE(algo.hasEulerianPath());
}

TEST_F(EulerianPathGTest, testLongPathDoesNotOverflowTheStack) {
    constexpr count n = 400000;
    Graph G(n);
    for (node u = 0; u + 1 < n; ++u) {
        G.addEdge(u, u + 1);
    }

    EulerianPath algo(G);
    algo.run();

    ASSERT_TRUE(algo.hasEulerianPath());
    EXPECT_EQ(algo.getPath().size(), n);
}

TEST_F(EulerianPathGTest, testRunTwice) {
    Graph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    EulerianPath algo(G);
    algo.run();
    const auto path = algo.getPath();
    algo.run();

    EXPECT_EQ(algo.getPath(), path);
}

TEST_F(EulerianPathGTest, testResultsBeforeRunThrow) {
    Graph G(2);
    EulerianPath algo(G);
    EXPECT_THROW(algo.getPath(), std::runtime_error);
    EXPECT_THROW(algo.hasEulerianPath(), std::runtime_error);
}

} /* namespace NetworKit */

// networkit-format
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/MinimumArborescence.hpp>

namespace NetworKit {

class MinimumArborescenceGTest : public testing::Test {
protected:
    static count numberOfReachable(const Graph &G, node root) {
        std::vector<bool> visited(G.upperNodeIdBound(), false);
        std::queue<node> queue;
        visited[root] = true;
        queue.push(root);
        count reached = 0;
        while (!queue.empty()) {
            const node u = queue.front();
            queue.pop();
            ++reached;
            G.forNeighborsOf(u, [&](node v) -> void {
                if (!visited[v]) {
                    visited[v] = true;
                    queue.push(v);
                }
            });
        }
        return reached;
    }

    static void expectValidArborescence(const Graph &G, const Graph &tree, node root) {
        EXPECT_EQ(tree.numberOfEdges(), G.numberOfNodes() - 1);
        EXPECT_EQ(tree.degreeIn(root), 0);
        EXPECT_EQ(numberOfReachable(tree, root), G.numberOfNodes());
        tree.forEdges([&](node u, node v, edgeweight w) -> void {
            EXPECT_TRUE(G.hasEdge(u, v));
            EXPECT_DOUBLE_EQ(G.weight(u, v), w);
        });
    }

    // Tries all assignments of an incoming edge to each node. Returns infinity if none of them is
    // an arborescence. Node ids must be consecutive.
    static edgeweight bruteForce(const Graph &G, node root) {
        const count n = G.numberOfNodes();
        edgeweight best = std::numeric_limits<edgeweight>::infinity();
        Graph tree(n, true, true);

        auto assign = [&](auto &&self, node u, edgeweight weight) -> void {
            if (u == n) {
                if (numberOfReachable(tree, root) == n && weight < best) {
                    best = weight;
                }
                return;
            }
            if (u == root) {
                self(self, u + 1, weight);
                return;
            }
            G.forInEdgesOf(u, [&](node, node v, edgeweight w) -> void {
                if (v != u) {
                    tree.addEdge(v, u, w);
                    self(self, u + 1, weight + w);
                    tree.removeEdge(v, u);
                }
            });
        };
        assign(assign, 0, 0.);
        return best;
    }
};

TEST_F(MinimumArborescenceGTest, testRequiresCycleContraction) {
    // The cheapest incoming edges of 1 and 2 form a cycle, one of them has to be given up.
    Graph G(3, true, true);
    G.addEdge(0, 1, 10);
    G.addEdge(0, 2, 11);
    G.addEdge(1, 2, 1);
    G.addEdge(2, 1, 1);

    MinimumArborescence algo(G, 0);
    algo.run();

    ASSERT_TRUE(algo.hasArborescence());
    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 11);
    expectValidArborescence(G, algo.getArborescence(), 0);
}

TEST_F(MinimumArborescenceGTest, testNestedCycles) {
    Graph G(5, true, true);
    G.addEdge(0, 1, 100);
    G.addEdge(1, 2, 1);
    G.addEdge(2, 1, 1);
    G.addEdge(2, 3, 2);
    G.addEdge(3, 4, 1);
    G.addEdge(4, 3, 1);
    G.addEdge(4, 1, 3);
    G.addEdge(0, 3, 50);

    MinimumArborescence algo(G, 0);
    algo.run();

    ASSERT_TRUE(algo.hasArborescence());
    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), bruteForce(G, 0));
    expectValidArborescence(G, algo.getArborescence(), 0);
}

TEST_F(MinimumArborescenceGTest, testUnreachableNode) {
    Graph G(3, false, true);
    G.addEdge(0, 1);
    G.addEdge(2, 1);

    MinimumArborescence algo(G, 0);
    algo.run();

    EXPECT_FALSE(algo.hasArborescence());
    EXPECT_EQ(algo.getArborescence().numberOfEdges(), 0);
    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 0);
}

TEST_F(MinimumArborescenceGTest, testReachabilityDependsOnTheRoot) {
    // All nodes have an incoming edge, but the cycle is not reachable from node 0.
    Graph G(4, false, true);
    G.addEdge(0, 3);
    G.addEdge(1, 2);
    G.addEdge(2, 1);

    MinimumArborescence algo(G, 0);
    algo.run();

    EXPECT_FALSE(algo.hasArborescence());
}

TEST_F(MinimumArborescenceGTest, testNegativeWeights) {
    Graph G(3, true, true);
    G.addEdge(0, 1, -5);
    G.addEdge(1, 2, -2);
    G.addEdge(0, 2, -1);

    MinimumArborescence algo(G, 0);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), -7);
}

TEST_F(MinimumArborescenceGTest, testSelfLoopsParallelEdgesAndEdgesIntoTheRoot) {
    Graph G(3, true, true);
    G.addEdge(0, 0, -100);
    G.addEdge(1, 1, -100);
    G.addEdge(1, 0, -100);
    G.addEdge(0, 1, 5);
    G.addEdge(0, 1, 3);
    G.addEdge(1, 2, 1);

    MinimumArborescence algo(G, 0);
    algo.run();

    ASSERT_TRUE(algo.hasArborescence());
    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 4);
    EXPECT_EQ(algo.getArborescence().numberOfEdges(), 2);
}

TEST_F(MinimumArborescenceGTest, testDeletedNodes) {
    Graph G(5, true, true);
    G.addEdge(1, 3, 2);
    G.addEdge(3, 4, 2);
    G.addEdge(1, 4, 5);
    G.addEdge(0, 3, 1);
    G.removeNode(0);
    G.removeNode(2);

    MinimumArborescence algo(G, 1);
    algo.run();

    ASSERT_TRUE(algo.hasArborescence());
    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 4);
    EXPECT_FALSE(algo.getArborescence().hasNode(0));
}

TEST_F(MinimumArborescenceGTest, testSingleNode) {
    Graph G(1, false, true);
    MinimumArborescence algo(G, 0);
    algo.run();
    EXPECT_TRUE(algo.hasArborescence());
    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 0);
}

TEST_F(MinimumArborescenceGTest, testMatchesBruteForce) {
    for (const count seed : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}) {
        Aux::Random::setSeed(seed, false);
        for (const double probability : {0.4, 0.7}) {
            constexpr count n = 7;
            Graph G(n, true, true);
            for (node u = 0; u < n; ++u) {
                for (node v = 0; v < n; ++v) {
                    if (u != v && v != 0 && Aux::Random::probability() < probability) {
                        G.addEdge(u, v, static_cast<edgeweight>(Aux::Random::integer(1, 9)));
                    }
                }
            }

            MinimumArborescence algo(G, 0);
            algo.run();

            const edgeweight expected = bruteForce(G, 0);
            if (expected == std::numeric_limits<edgeweight>::infinity()) {
                EXPECT_FALSE(algo.hasArborescence());
            } else {
                ASSERT_TRUE(algo.hasArborescence());
                EXPECT_DOUBLE_EQ(algo.getTotalWeight(), expected);
                expectValidArborescence(G, algo.getArborescence(), 0);
            }
        }
    }
}

TEST_F(MinimumArborescenceGTest, testRunTwice) {
    Graph G(3, true, true);
    G.addEdge(0, 1, 2);
    G.addEdge(1, 2, 3);

    MinimumArborescence algo(G, 0);
    algo.run();
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getTotalWeight(), 5);
    EXPECT_EQ(algo.getArborescence().numberOfEdges(), 2);
}

TEST_F(MinimumArborescenceGTest, testInvalidInput) {
    EXPECT_THROW(MinimumArborescence(Graph(2), 0), std::runtime_error);
    EXPECT_THROW(MinimumArborescence(Graph(2, false, true), 5), std::runtime_error);
}

TEST_F(MinimumArborescenceGTest, testResultsBeforeRunThrow) {
    Graph G(2, false, true);
    MinimumArborescence algo(G, 0);
    EXPECT_THROW(algo.getArborescence(), std::runtime_error);
    EXPECT_THROW(algo.getTotalWeight(), std::runtime_error);
    EXPECT_THROW(algo.hasArborescence(), std::runtime_error);
}

} /* namespace NetworKit */

// networkit-format
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/flow/GomoryHuTree.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class GomoryHuTreeGTest : public testing::Test {
protected:
    static Graph randomWeightedGraph(count n, double probability) {
        Graph G = GraphTools::toWeighted(ErdosRenyiGenerator(n, probability).generate());
        G.forEdges([&](node u, node v) -> void {
            G.setWeight(u, v, static_cast<edgeweight>(Aux::Random::integer(1, 9)));
        });
        return G;
    }

    // Total weight of the edges between the nodes in the bit mask and the others.
    static edgeweight cutWeight(const Graph &G, count mask) {
        edgeweight weight = 0;
        G.forEdges([&](node u, node v, edgeweight w) -> void {
            if (((mask >> u) & 1) != ((mask >> v) & 1)) {
                weight += w;
            }
        });
        return weight;
    }

    // Minimum cut between s and t, enumerating all cuts of a graph with consecutive node ids.
    static edgeweight bruteForceMinCut(const Graph &G, node s, node t) {
        edgeweight best = std::numeric_limits<edgeweight>::infinity();
        for (count mask = 0; mask < (count{1} << G.numberOfNodes()); ++mask) {
            if (((mask >> s) & 1) && !((mask >> t) & 1)) {
                best = std::min(best, cutWeight(G, mask));
            }
        }
        return best;
    }

    // Bit mask of the nodes reachable from start in the tree without using the edge (start, skip).
    static count treeSide(const Graph &tree, node start, node skip) {
        count mask = count{1} << start;
        std::vector<node> stack{start};
        while (!stack.empty()) {
            const node u = stack.back();
            stack.pop_back();
            tree.forNeighborsOf(u, [&](node v) -> void {
                if (!((mask >> v) & 1) && !(u == start && v == skip)) {
                    mask |= count{1} << v;
                    stack.push_back(v);
                }
            });
        }
        return mask;
    }
};

TEST_F(GomoryHuTreeGTest, testPath) {
    Graph G(4, true);
    G.addEdge(0, 1, 3);
    G.addEdge(1, 2, 1);
    G.addEdge(2, 3, 2);

    GomoryHuTree algo(G);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getMinCutValue(0, 1), 3);
    EXPECT_DOUBLE_EQ(algo.getMinCutValue(0, 3), 1);
    EXPECT_DOUBLE_EQ(algo.getMinCutValue(3, 2), 2);
    EXPECT_EQ(algo.getTree().numberOfEdges(), 3);
}

TEST_F(GomoryHuTreeGTest, testMatchesBruteForce) {
    for (const count seed : {1, 2, 3, 4, 5, 6}) {
        Aux::Random::setSeed(seed, false);
        for (const double probability : {0.3, 0.6, 0.9}) {
            const Graph G = randomWeightedGraph(9, probability);

            GomoryHuTree algo(G);
            algo.run();

            for (node u = 0; u < 9; ++u) {
                for (node v = u + 1; v < 9; ++v) {
                    EXPECT_DOUBLE_EQ(algo.getMinCutValue(u, v), bruteForceMinCut(G, u, v));
                }
            }
        }
    }
}

TEST_F(GomoryHuTreeGTest, testTreeEdgesInduceMinimumCuts) {
    for (const count seed : {1, 2, 3, 4, 5, 6}) {
        Aux::Random::setSeed(seed, false);
        const Graph G = randomWeightedGraph(12, 0.5);

        GomoryHuTree algo(G);
        algo.run();

        const Graph tree = algo.getTree();
        EXPECT_EQ(tree.numberOfEdges(), 11);
        tree.forEdges([&](node u, node v, edgeweight w) -> void {
            EXPECT_DOUBLE_EQ(cutWeight(G, treeSide(tree, u, v)), w);
        });
    }
}

TEST_F(GomoryHuTreeGTest, testDisconnectedGraph) {
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(2, 3);

    GomoryHuTree algo(G);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getMinCutValue(0, 1), 1);
    EXPECT_DOUBLE_EQ(algo.getMinCutValue(0, 2), 0);
    EXPECT_EQ(algo.getTree().numberOfEdges(), 3);
}

TEST_F(GomoryHuTreeGTest, testSelfLoopsAndParallelEdges) {
    Graph G(3, true);
    G.addEdge(0, 0, 100);
    G.addEdge(0, 1, 1);
    G.addEdge(0, 1, 2);
    G.addEdge(1, 2, 4);

    GomoryHuTree algo(G);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getMinCutValue(0, 1), 3);
    EXPECT_DOUBLE_EQ(algo.getMinCutValue(1, 2), 4);
}

TEST_F(GomoryHuTreeGTest, testDeletedNodes) {
    Graph G(5);
    G.addEdge(1, 3);
    G.addEdge(3, 4);
    G.addEdge(4, 1);
    G.addEdge(0, 3);
    G.removeNode(0);

    GomoryHuTree algo(G);
    algo.run();

    const Graph tree = algo.getTree();
    EXPECT_EQ(tree.numberOfNodes(), 4);
    EXPECT_FALSE(tree.hasNode(0));
    EXPECT_DOUBLE_EQ(algo.getMinCutValue(1, 4), 2);
    EXPECT_DOUBLE_EQ(algo.getMinCutValue(1, 2), 0);
    EXPECT_THROW(algo.getMinCutValue(0, 1), std::runtime_error);
}

TEST_F(GomoryHuTreeGTest, testEmptyAndSingleNodeGraphs) {
    Graph empty;
    GomoryHuTree onEmpty(empty);
    onEmpty.run();
    EXPECT_EQ(onEmpty.getTree().numberOfNodes(), 0);

    Graph single(1);
    GomoryHuTree onSingle(single);
    onSingle.run();
    EXPECT_EQ(onSingle.getTree().numberOfEdges(), 0);
}

TEST_F(GomoryHuTreeGTest, testRunTwice) {
    Graph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    GomoryHuTree algo(G);
    algo.run();
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getMinCutValue(0, 2), 1);
    EXPECT_EQ(algo.getTree().numberOfEdges(), 2);
}

TEST_F(GomoryHuTreeGTest, testInvalidInput) {
    EXPECT_THROW(GomoryHuTree{Graph(3, false, true)}, std::runtime_error);

    Graph negative(2, true);
    negative.addEdge(0, 1, -1);
    EXPECT_THROW(GomoryHuTree{negative}, std::runtime_error);

    Graph G(2);
    GomoryHuTree algo(G);
    EXPECT_THROW(algo.getTree(), std::runtime_error);
    EXPECT_THROW(algo.getMinCutValue(0, 1), std::runtime_error);
    algo.run();
    EXPECT_THROW(algo.getMinCutValue(0, 0), std::runtime_error);
}

} // namespace NetworKit

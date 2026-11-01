// networkit-format
#include <set>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/community/Modularity.hpp>
#include <networkit/community/Walktrap.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class WalktrapGTest : public testing::Test {
protected:
    static void addClique(Graph &G, node first, count size, edgeweight weight = 1) {
        for (node u = first; u < first + size; ++u) {
            for (node v = u + 1; v < first + size; ++v) {
                G.addEdge(u, v, weight);
            }
        }
    }

    static Graph ringOfCliques(count cliques, count size) {
        Graph G(cliques * size);
        for (index i = 0; i < cliques; ++i) {
            addClique(G, i * size, size);
            G.addEdge(i * size, ((i + 1) % cliques) * size + 1);
        }
        return G;
    }

    static void expectBlocks(const Partition &zeta, const std::vector<std::vector<node>> &blocks) {
        std::set<index> ids;
        for (const auto &block : blocks) {
            ids.insert(zeta[block.front()]);
            for (const node u : block) {
                EXPECT_EQ(zeta[u], zeta[block.front()]);
            }
        }
        EXPECT_EQ(ids.size(), blocks.size());
    }
};

TEST_F(WalktrapGTest, testRingOfCliques) {
    constexpr count cliques = 6;
    constexpr count cliqueSize = 5;
    const Graph G = ringOfCliques(cliques, cliqueSize);

    Walktrap algo(G);
    algo.run();

    const Partition &zeta = algo.getPartition();
    EXPECT_EQ(zeta.numberOfSubsets(), cliques);
    for (node u = 0; u < G.numberOfNodes(); ++u) {
        EXPECT_EQ(zeta[u], zeta[u - u % cliqueSize]);
    }
}

TEST_F(WalktrapGTest, testWalkLengthsFindTheSameClearStructure) {
    constexpr count cliques = 4;
    constexpr count cliqueSize = 6;
    const Graph G = ringOfCliques(cliques, cliqueSize);

    for (const count walkLength : {2, 3, 4, 8}) {
        Walktrap algo(G, walkLength);
        algo.run();
        EXPECT_EQ(algo.getPartition().numberOfSubsets(), cliques);
    }
}

TEST_F(WalktrapGTest, testWeightsDecideTheCommunities) {
    // All pairs are adjacent, so only the weights tell the two groups apart.
    Graph G(6, true);
    addClique(G, 0, 6, 0.1);
    for (const node offset : {0, 3}) {
        G.setWeight(offset, offset + 1, 10);
        G.setWeight(offset, offset + 2, 10);
        G.setWeight(offset + 1, offset + 2, 10);
    }

    Walktrap algo(G);
    algo.run();

    expectBlocks(algo.getPartition(), {{0, 1, 2}, {3, 4, 5}});
}

TEST_F(WalktrapGTest, testScalingTheWeightsChangesNothing) {
    Graph G(12, true);
    addClique(G, 0, 6, 2);
    addClique(G, 6, 6, 2);
    G.addEdge(0, 6, 0.5);
    Graph scaled = G;
    scaled.forEdges(
        [&](node u, node v, edgeweight weight) -> void { scaled.setWeight(u, v, weight * 7); });

    Walktrap algo(G), scaledAlgo(scaled);
    algo.run();
    scaledAlgo.run();

    EXPECT_EQ(algo.getPartition().getVector(), scaledAlgo.getPartition().getVector());
}

TEST_F(WalktrapGTest, testComponentsAreNeverMerged) {
    // Two edges: merging them would lower the modularity, but it is impossible regardless.
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(2, 3);

    Walktrap algo(G);
    algo.run();

    expectBlocks(algo.getPartition(), {{0, 1}, {2, 3}});
}

TEST_F(WalktrapGTest, testIsolatedNodesStaySingletons) {
    Graph G(7);
    addClique(G, 0, 4);

    Walktrap algo(G);
    algo.run();

    expectBlocks(algo.getPartition(), {{0, 1, 2, 3}, {4}, {5}, {6}});
}

TEST_F(WalktrapGTest, testEdgelessAndEmptyGraphs) {
    Graph edgeless(3), empty;

    Walktrap algo(edgeless), emptyAlgo(empty);
    algo.run();
    emptyAlgo.run();

    EXPECT_EQ(algo.getPartition().numberOfSubsets(), 3);
    EXPECT_EQ(emptyAlgo.getPartition().numberOfSubsets(), 0);
}

TEST_F(WalktrapGTest, testZeroWeightEdgesDoNotConnect) {
    Graph G(4, true);
    G.addEdge(0, 1, 1);
    G.addEdge(2, 3, 1);
    G.addEdge(1, 2, 0);

    Walktrap algo(G);
    algo.run();

    expectBlocks(algo.getPartition(), {{0, 1}, {2, 3}});
}

TEST_F(WalktrapGTest, testSelfLoops) {
    Graph G(10);
    addClique(G, 0, 5);
    addClique(G, 5, 5);
    G.addEdge(0, 5);
    G.addEdge(0, 0);
    G.addEdge(7, 7);

    Walktrap algo(G);
    algo.run();

    expectBlocks(algo.getPartition(), {{0, 1, 2, 3, 4}, {5, 6, 7, 8, 9}});
}

TEST_F(WalktrapGTest, testDeletedNodes) {
    Graph G(12);
    addClique(G, 0, 6);
    addClique(G, 6, 6);
    G.addEdge(0, 6);
    G.removeNode(3);
    G.removeNode(9);

    Walktrap algo(G);
    algo.run();

    const Partition &zeta = algo.getPartition();
    expectBlocks(zeta, {{0, 1, 2, 4, 5}, {6, 7, 8, 10, 11}});
    EXPECT_EQ(zeta.numberOfSubsets(), 2);
    EXPECT_EQ(zeta[3], none);
    EXPECT_EQ(zeta[9], none);
}

TEST_F(WalktrapGTest, testKarateClub) {
    // Estimating merge costs, as igraph does, gives 5 communities with a modularity of 0.3532.
    const Graph G = METISGraphReader{}.read("input/karate.graph");

    Walktrap algo(G);
    algo.run();

    EXPECT_EQ(algo.getPartition().numberOfSubsets(), 4);
    EXPECT_NEAR(Modularity{}.getQuality(algo.getPartition(), G), 0.4198, 1e-4);
}

TEST_F(WalktrapGTest, testRunTwice) {
    const Graph G = ringOfCliques(3, 5);

    Walktrap algo(G);
    algo.run();
    const std::vector<index> first = algo.getPartition().getVector();
    algo.run();

    EXPECT_EQ(algo.getPartition().getVector(), first);
}

TEST_F(WalktrapGTest, testInvalidInput) {
    EXPECT_THROW(Walktrap(Graph(2, false, true)), std::runtime_error);
    EXPECT_THROW(Walktrap(Graph(2), 0), std::runtime_error);

    Graph negative(2, true);
    negative.addEdge(0, 1, -1);
    EXPECT_THROW(Walktrap{negative}, std::runtime_error);
}

TEST_F(WalktrapGTest, testPartitionBeforeRunThrows) {
    Graph G(2);
    Walktrap algo(G);
    EXPECT_THROW(algo.getPartition(), std::runtime_error);
}

} /* namespace NetworKit */

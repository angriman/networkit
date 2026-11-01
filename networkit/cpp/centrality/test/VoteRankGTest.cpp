// networkit-format
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/centrality/VoteRank.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class VoteRankGTest : public testing::Test {};

TEST_F(VoteRankGTest, testToyGraph) {
    Graph G(10);
    G.addEdge(7, 8);
    G.addEdge(7, 5);
    G.addEdge(7, 9);
    G.addEdge(5, 0);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(0, 3);
    G.addEdge(0, 4);
    G.addEdge(1, 6);
    G.addEdge(2, 6);
    G.addEdge(3, 6);
    G.addEdge(4, 6);

    VoteRank algo(G, 3);
    algo.run();

    EXPECT_EQ(algo.getSpreaders(), std::vector<node>({0, 7, 6}));
}

TEST_F(VoteRankGTest, testStopsWhenNoNodeGetsVotes) {
    Graph G(5);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(0, 3);
    G.addEdge(1, 4);

    VoteRank algo(G);
    algo.run();

    EXPECT_EQ(algo.getSpreaders(), std::vector<node>({0, 1}));
}

TEST_F(VoteRankGTest, testDirected) {
    // Node 0 has the most out-neighbors, hence it receives the most votes.
    Graph G(4, false, true);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(0, 3);
    G.addEdge(1, 2);

    VoteRank algo(G, 1);
    algo.run();

    EXPECT_EQ(algo.getSpreaders(), std::vector<node>({0}));
}

TEST_F(VoteRankGTest, testDeletedNodes) {
    Graph G(5);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(3, 4);
    G.removeNode(4);

    VoteRank algo(G, 1);
    algo.run();

    EXPECT_EQ(algo.getSpreaders(), std::vector<node>({1}));
}

TEST_F(VoteRankGTest, testGraphWithoutEdges) {
    Graph G(3);
    VoteRank algo(G);
    algo.run();
    EXPECT_TRUE(algo.getSpreaders().empty());
}

TEST_F(VoteRankGTest, testEmptyGraph) {
    Graph G;
    VoteRank algo(G);
    algo.run();
    EXPECT_TRUE(algo.getSpreaders().empty());
}

TEST_F(VoteRankGTest, testZeroSpreaders) {
    Graph G(2);
    G.addEdge(0, 1);
    VoteRank algo(G, 0);
    algo.run();
    EXPECT_TRUE(algo.getSpreaders().empty());
}

TEST_F(VoteRankGTest, testRunTwice) {
    Graph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    VoteRank algo(G);
    algo.run();
    const auto spreaders = algo.getSpreaders();
    algo.run();

    EXPECT_EQ(algo.getSpreaders(), spreaders);
}

TEST_F(VoteRankGTest, testResultsBeforeRunThrow) {
    Graph G(2);
    VoteRank algo(G);
    EXPECT_THROW(algo.getSpreaders(), std::runtime_error);
}

} /* namespace NetworKit */

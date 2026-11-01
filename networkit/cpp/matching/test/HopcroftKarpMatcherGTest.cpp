// networkit-format
#include <functional>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/matching/HopcroftKarpMatcher.hpp>

namespace NetworKit {

class HopcroftKarpMatcherGTest : public testing::Test {
protected:
    // Random bipartite graph, the nodes [0, numLeft) form one side.
    static Graph randomBipartite(count numLeft, count numRight, double probability) {
        Graph G(numLeft + numRight);
        for (node u = 0; u < numLeft; ++u) {
            for (node v = numLeft; v < numLeft + numRight; ++v) {
                if (Aux::Random::probability() < probability) {
                    G.addEdge(u, v);
                }
            }
        }
        return G;
    }

    // Size of a maximum matching computed with Kuhn's algorithm, as an independent reference.
    static count kuhnMatchingSize(const Graph &G, count numLeft) {
        std::vector<node> mate(G.upperNodeIdBound(), none);
        std::vector<bool> visited;

        std::function<bool(node)> augment = [&](node u) -> bool {
            for (const node v : G.neighborRange(u)) {
                if (visited[v]) {
                    continue;
                }
                visited[v] = true;
                if (mate[v] == none || augment(mate[v])) {
                    mate[v] = u;
                    return true;
                }
            }
            return false;
        };

        count size = 0;
        for (node u = 0; u < numLeft; ++u) {
            visited.assign(G.upperNodeIdBound(), false);
            if (augment(u)) {
                ++size;
            }
        }
        return size;
    }
};

TEST_F(HopcroftKarpMatcherGTest, testPerfectMatchingRequiresAugmenting) {
    // A greedy matching picking (1, 4) and (0, 3) first gets stuck.
    Graph G(6);
    G.addEdge(0, 3);
    G.addEdge(1, 3);
    G.addEdge(1, 4);
    G.addEdge(2, 4);
    G.addEdge(2, 5);

    HopcroftKarpMatcher matcher(G);
    matcher.run();

    const auto matching = matcher.getMatching();
    EXPECT_TRUE(matching.isProper(G));
    EXPECT_EQ(matching.size(G), 3);
}

TEST_F(HopcroftKarpMatcherGTest, testStar) {
    Graph G(5);
    for (node u = 1; u < 5; ++u) {
        G.addEdge(0, u);
    }

    HopcroftKarpMatcher matcher(G);
    matcher.run();

    EXPECT_EQ(matcher.getMatching().size(G), 1);
}

TEST_F(HopcroftKarpMatcherGTest, testDisconnectedGraphWithIsolatedNodes) {
    Graph G(7);
    G.addEdge(0, 1);
    G.addEdge(3, 4);
    G.addEdge(4, 5);

    HopcroftKarpMatcher matcher(G);
    matcher.run();

    const auto matching = matcher.getMatching();
    EXPECT_TRUE(matching.isProper(G));
    EXPECT_EQ(matching.size(G), 2);
    EXPECT_FALSE(matching.isMatched(2));
    EXPECT_FALSE(matching.isMatched(6));
}

TEST_F(HopcroftKarpMatcherGTest, testDeletedNodes) {
    Graph G(5);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.removeNode(2);

    HopcroftKarpMatcher matcher(G);
    matcher.run();

    const auto matching = matcher.getMatching();
    EXPECT_TRUE(matching.isProper(G));
    EXPECT_EQ(matching.size(G), 2);
}

TEST_F(HopcroftKarpMatcherGTest, testEmptyGraph) {
    Graph G;
    HopcroftKarpMatcher matcher(G);
    matcher.run();
    EXPECT_EQ(matcher.getMatching().size(G), 0);
}

TEST_F(HopcroftKarpMatcherGTest, testParallelEdgesAndWeightsAreIgnored) {
    Graph G(4, true);
    G.addEdge(0, 1, 100.);
    G.addEdge(0, 1, 100.);
    G.addEdge(1, 2, 1.);
    G.addEdge(2, 3, 1.);

    HopcroftKarpMatcher matcher(G);
    matcher.run();

    EXPECT_EQ(matcher.getMatching().size(G), 2);
}

TEST_F(HopcroftKarpMatcherGTest, testLongAugmentingPath) {
    // Worst case for a recursive search. The path n - 2, n - 3, ..., 0, n - 1 has its edges added
    // such that the first phase matches each even node u to u + 1. Only the ends of the path stay
    // free, connected by a single augmenting path through all nodes.
    constexpr count n = 400000;
    Graph G(n);
    for (node u = 0; u + 2 < n; u += 2) {
        G.addEdge(u, u + 1);
    }
    for (node u = 1; u + 2 < n; u += 2) {
        G.addEdge(u, u + 1);
    }
    G.addEdge(0, n - 1);

    HopcroftKarpMatcher matcher(G);
    matcher.run();

    const auto matching = matcher.getMatching();
    EXPECT_EQ(matching.size(G), n / 2);
    EXPECT_EQ(matching.mate(n - 1), 0);
}

TEST_F(HopcroftKarpMatcherGTest, testMatchesKuhnOnRandomGraphs) {
    for (const count seed : {1, 2, 3, 4, 5}) {
        Aux::Random::setSeed(seed, false);
        for (const double probability : {0.02, 0.05, 0.1, 0.3}) {
            constexpr count numLeft = 60, numRight = 80;
            const Graph G = randomBipartite(numLeft, numRight, probability);

            HopcroftKarpMatcher matcher(G);
            matcher.run();

            const auto matching = matcher.getMatching();
            EXPECT_TRUE(matching.isProper(G));
            EXPECT_EQ(matching.size(G), kuhnMatchingSize(G, numLeft));
        }
    }
}

TEST_F(HopcroftKarpMatcherGTest, testRunTwice) {
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);

    HopcroftKarpMatcher matcher(G);
    matcher.run();
    matcher.run();

    const auto matching = matcher.getMatching();
    EXPECT_TRUE(matching.isProper(G));
    EXPECT_EQ(matching.size(G), 2);
}

TEST_F(HopcroftKarpMatcherGTest, testInvalidInput) {
    Graph directed(2, false, true);
    EXPECT_THROW(HopcroftKarpMatcher{directed}, std::invalid_argument);

    Graph triangle(3);
    triangle.addEdge(0, 1);
    triangle.addEdge(1, 2);
    triangle.addEdge(2, 0);
    EXPECT_THROW(HopcroftKarpMatcher{triangle}, std::invalid_argument);

    Graph selfLoop(2);
    selfLoop.addEdge(0, 0);
    EXPECT_THROW(HopcroftKarpMatcher{selfLoop}, std::invalid_argument);
}

} /* namespace NetworKit */

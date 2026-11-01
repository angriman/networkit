// networkit-format
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/matching/BlossomMatcher.hpp>

namespace NetworKit {

class BlossomMatcherGTest : public testing::Test {
protected:
    static Matching matchingOf(const Graph &G) {
        BlossomMatcher matcher(G);
        matcher.run();
        const Matching matching = matcher.getMatching();
        EXPECT_TRUE(matching.isProper(G));
        return matching;
    }

    // Exhaustive search over the subsets of nodes. Node ids must be consecutive and few.
    static edgeweight optimalWeight(const Graph &G) {
        const count n = G.numberOfNodes();
        std::vector<edgeweight> best(uint64_t{1} << n, 0);
        for (uint64_t subset = 1; subset < best.size(); ++subset) {
            node u = 0;
            while (((subset >> u) & 1) == 0) {
                ++u;
            }
            const uint64_t rest = subset & ~(uint64_t{1} << u);
            best[subset] = best[rest];
            G.forEdgesOf(u, [&](node, node v, edgeweight weight) -> void {
                if (v != u && ((rest >> v) & 1) != 0) {
                    best[subset] =
                        std::max(best[subset], weight + best[rest & ~(uint64_t{1} << v)]);
                }
            });
        }
        return best.back();
    }

    static Graph randomGraph(count n, double probability, count maxWeight) {
        Graph G(n, true);
        for (node u = 0; u < n; ++u) {
            for (node v = u + 1; v < n; ++v) {
                if (Aux::Random::probability() < probability) {
                    G.addEdge(u, v, static_cast<edgeweight>(Aux::Random::integer(1, maxWeight)));
                }
            }
        }
        return G;
    }
};

TEST_F(BlossomMatcherGTest, testOddCycleWithPendantEdges) {
    Graph G(5, true);
    G.addEdge(0, 1, 6);
    G.addEdge(1, 2, 5);
    G.addEdge(2, 0, 5);
    G.addEdge(1, 3, 4);
    G.addEdge(2, 4, 4);

    EXPECT_DOUBLE_EQ(matchingOf(G).weight(G), 10);
}

TEST_F(BlossomMatcherGTest, testPrefersWeightOverCardinality) {
    Graph G(4, true);
    G.addEdge(0, 1, 1);
    G.addEdge(1, 2, 10);
    G.addEdge(2, 3, 1);

    const Matching matching = matchingOf(G);

    EXPECT_DOUBLE_EQ(matching.weight(G), 10);
    EXPECT_EQ(matching.mate(1), 2);
}

TEST_F(BlossomMatcherGTest, testUnweightedGraphsGetMaximumCardinality) {
    // Two triangles joined by an edge, which a perfect matching has to use.
    Graph G(6);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 0);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 3);
    G.addEdge(2, 3);

    const Matching matching = matchingOf(G);

    EXPECT_EQ(matching.size(G), 3);
    EXPECT_EQ(matching.mate(2), 3);
}

TEST_F(BlossomMatcherGTest, testMatchesExhaustiveSearch) {
    for (const count seed : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}) {
        Aux::Random::setSeed(seed, false);
        for (const double probability : {0.2, 0.5, 0.9}) {
            for (const count maxWeight : {1, 3, 20}) {
                const Graph G = randomGraph(12, probability, maxWeight);
                EXPECT_DOUBLE_EQ(matchingOf(G).weight(G), optimalWeight(G));
            }
        }
    }
}

TEST_F(BlossomMatcherGTest, testFractionalWeights) {
    Aux::Random::setSeed(11, false);
    for (count i = 0; i < 20; ++i) {
        Graph G(10, true);
        for (node u = 0; u < 10; ++u) {
            for (node v = u + 1; v < 10; ++v) {
                if (Aux::Random::probability() < 0.6) {
                    G.addEdge(u, v, Aux::Random::probability());
                }
            }
        }
        EXPECT_NEAR(matchingOf(G).weight(G), optimalWeight(G), 1e-9);
    }
}

TEST_F(BlossomMatcherGTest, testNonPositiveWeightsAreNeverMatched) {
    Graph G(4, true);
    G.addEdge(0, 1, -3);
    G.addEdge(1, 2, 0);
    G.addEdge(2, 3, 2);

    const Matching matching = matchingOf(G);

    EXPECT_EQ(matching.size(G), 1);
    EXPECT_EQ(matching.mate(0), none);
    EXPECT_EQ(matching.mate(2), 3);
}

TEST_F(BlossomMatcherGTest, testSelfLoopsAndParallelEdges) {
    Graph G(3, true);
    G.addEdge(0, 0, 100);
    G.addEdge(1, 2, 1);
    G.addEdge(1, 2, 5);
    G.addEdge(0, 1, 2);

    const Matching matching = matchingOf(G);

    EXPECT_EQ(matching.mate(1), 2);
    EXPECT_EQ(matching.mate(0), none);
}

TEST_F(BlossomMatcherGTest, testDeletedNodes) {
    Graph G(6, true);
    G.addEdge(1, 3, 2);
    G.addEdge(3, 5, 3);
    G.addEdge(0, 1, 9);
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);

    const Matching matching = matchingOf(G);

    EXPECT_EQ(matching.mate(3), 5);
    EXPECT_EQ(matching.mate(1), none);
}

TEST_F(BlossomMatcherGTest, testEmptyAndEdgelessGraphs) {
    EXPECT_EQ(matchingOf(Graph()).size(Graph()), 0);

    const Graph edgeless(4);
    EXPECT_EQ(matchingOf(edgeless).size(edgeless), 0);
}

TEST_F(BlossomMatcherGTest, testRunTwice) {
    Graph G(2, true);
    G.addEdge(0, 1, 2);

    BlossomMatcher matcher(G);
    matcher.run();
    matcher.run();

    EXPECT_EQ(matcher.getMatching().mate(0), 1);
}

TEST_F(BlossomMatcherGTest, testInvalidInput) {
    EXPECT_THROW(BlossomMatcher{Graph(2, false, true)}, std::invalid_argument);
}

TEST_F(BlossomMatcherGTest, testResultBeforeRunThrows) {
    Graph G(2);
    BlossomMatcher matcher(G);
    EXPECT_THROW(matcher.getMatching(), std::runtime_error);
}

} /* namespace NetworKit */

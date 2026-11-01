// networkit-format
#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/global/GraphletCensus.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

using Graphlet = GraphletCensus::Graphlet;

class GraphletCensusGTest : public testing::Test {
protected:
    static std::vector<count> expectedOnly(Graphlet graphlet, count occurrences = 1) {
        std::vector<count> expected(GraphletCensus::numberOfGraphlets, 0);
        expected[static_cast<index>(graphlet)] = occurrences;
        return expected;
    }

    static std::vector<count> censusOf(const Graph &G) {
        GraphletCensus census(G);
        census.run();
        return census.getCounts();
    }

    static std::vector<count> bruteForce(const Graph &G) {
        std::vector<count> census(GraphletCensus::numberOfGraphlets, 0);
        const auto increment = [&](Graphlet graphlet) -> void {
            ++census[static_cast<index>(graphlet)];
        };
        const count n = G.numberOfNodes();

        for (node a = 0; a < n; ++a) {
            for (node b = a + 1; b < n; ++b) {
                for (node c = b + 1; c < n; ++c) {
                    const count edges = G.hasEdge(a, b) + G.hasEdge(a, c) + G.hasEdge(b, c);
                    if (edges == 2) {
                        increment(Graphlet::Wedge);
                    } else if (edges == 3) {
                        increment(Graphlet::Triangle);
                    }

                    for (node d = c + 1; d < n; ++d) {
                        const std::array<node, 4> nodes{a, b, c, d};
                        std::array<count, 4> degrees{};
                        for (index i = 0; i < 4; ++i) {
                            for (index j = 0; j < 4; ++j) {
                                degrees[i] += i != j && G.hasEdge(nodes[i], nodes[j]);
                            }
                        }
                        std::sort(degrees.begin(), degrees.end());
                        if (degrees == std::array<count, 4>{1, 1, 2, 2}) {
                            increment(Graphlet::Path);
                        } else if (degrees == std::array<count, 4>{1, 1, 1, 3}) {
                            increment(Graphlet::Star);
                        } else if (degrees == std::array<count, 4>{2, 2, 2, 2}) {
                            increment(Graphlet::Cycle);
                        } else if (degrees == std::array<count, 4>{1, 2, 2, 3}) {
                            increment(Graphlet::Paw);
                        } else if (degrees == std::array<count, 4>{2, 2, 3, 3}) {
                            increment(Graphlet::Diamond);
                        } else if (degrees == std::array<count, 4>{3, 3, 3, 3}) {
                            increment(Graphlet::Clique);
                        }
                    }
                }
            }
        }
        return census;
    }
};

TEST_F(GraphletCensusGTest, testGraphletsOnFourNodes) {
    Graph path(4);
    path.addEdge(0, 1);
    path.addEdge(1, 2);
    path.addEdge(2, 3);
    std::vector<count> expected = expectedOnly(Graphlet::Path);
    expected[static_cast<index>(Graphlet::Wedge)] = 2;
    EXPECT_EQ(censusOf(path), expected);

    Graph star(4);
    star.addEdge(0, 1);
    star.addEdge(0, 2);
    star.addEdge(0, 3);
    expected = expectedOnly(Graphlet::Star);
    expected[static_cast<index>(Graphlet::Wedge)] = 3;
    EXPECT_EQ(censusOf(star), expected);

    Graph cycle(4);
    cycle.addEdge(0, 1);
    cycle.addEdge(1, 2);
    cycle.addEdge(2, 3);
    cycle.addEdge(3, 0);
    expected = expectedOnly(Graphlet::Cycle);
    expected[static_cast<index>(Graphlet::Wedge)] = 4;
    EXPECT_EQ(censusOf(cycle), expected);

    Graph paw(4);
    paw.addEdge(0, 1);
    paw.addEdge(1, 2);
    paw.addEdge(2, 0);
    paw.addEdge(2, 3);
    expected = expectedOnly(Graphlet::Paw);
    expected[static_cast<index>(Graphlet::Wedge)] = 2;
    expected[static_cast<index>(Graphlet::Triangle)] = 1;
    EXPECT_EQ(censusOf(paw), expected);

    Graph diamond = cycle;
    diamond.addEdge(0, 2);
    expected = expectedOnly(Graphlet::Diamond);
    expected[static_cast<index>(Graphlet::Wedge)] = 2;
    expected[static_cast<index>(Graphlet::Triangle)] = 2;
    EXPECT_EQ(censusOf(diamond), expected);
}

TEST_F(GraphletCensusGTest, testCompleteGraph) {
    Graph G(6);
    for (node u = 0; u < 6; ++u) {
        for (node v = u + 1; v < 6; ++v) {
            G.addEdge(u, v);
        }
    }

    std::vector<count> expected = expectedOnly(Graphlet::Clique, 15);
    expected[static_cast<index>(Graphlet::Triangle)] = 20;
    EXPECT_EQ(censusOf(G), expected);
}

TEST_F(GraphletCensusGTest, testMatchesBruteForce) {
    for (const count seed : {1, 2, 3, 4, 5}) {
        Aux::Random::setSeed(seed, false);
        for (const double probability : {0.05, 0.2, 0.5, 0.9}) {
            constexpr count n = 20;
            Graph G(n);
            for (node u = 0; u < n; ++u) {
                for (node v = u + 1; v < n; ++v) {
                    if (Aux::Random::probability() < probability) {
                        G.addEdge(u, v);
                    }
                }
            }

            EXPECT_EQ(censusOf(G), bruteForce(G));
        }
    }
}

TEST_F(GraphletCensusGTest, testSelfLoopsParallelEdgesAndWeightsAreIgnored) {
    Graph G(3, true);
    G.addEdge(0, 1, 5);
    G.addEdge(0, 1, 2);
    G.addEdge(1, 1, 3);
    G.addEdge(1, 2, 0.5);

    EXPECT_EQ(censusOf(G), expectedOnly(Graphlet::Wedge));
}

TEST_F(GraphletCensusGTest, testDeletedNodes) {
    Graph G(6);
    G.addEdge(1, 3);
    G.addEdge(3, 5);
    G.addEdge(5, 1);
    G.addEdge(0, 1);
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);

    EXPECT_EQ(censusOf(G), expectedOnly(Graphlet::Triangle));
}

TEST_F(GraphletCensusGTest, testEmptyAndEdgelessGraphs) {
    const std::vector<count> nothing(GraphletCensus::numberOfGraphlets, 0);
    EXPECT_EQ(censusOf(Graph()), nothing);
    EXPECT_EQ(censusOf(Graph(5)), nothing);
}

TEST_F(GraphletCensusGTest, testGetCount) {
    Graph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 0);

    GraphletCensus census(G);
    census.run();
    census.run();

    EXPECT_EQ(census.getCount(Graphlet::Triangle), 1);
    EXPECT_EQ(census.getCount(Graphlet::Wedge), 0);
}

TEST_F(GraphletCensusGTest, testInvalidInput) {
    EXPECT_THROW(GraphletCensus{Graph(3, false, true)}, std::runtime_error);
}

TEST_F(GraphletCensusGTest, testResultsBeforeRunThrow) {
    Graph G(3);
    GraphletCensus census(G);
    EXPECT_THROW(census.getCount(Graphlet::Clique), std::runtime_error);
    EXPECT_THROW(census.getCounts(), std::runtime_error);
}

} /* namespace NetworKit */

// networkit-format
#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/global/TriadCensus.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

using Triad = TriadCensus::Triad;
using Edges = std::vector<std::pair<node, node>>;

class TriadCensusGTest : public testing::Test {
protected:
    // The definitions of the triads, independent of the code table of the algorithm.
    static const std::array<Edges, TriadCensus::numberOfTriads> &definitions() {
        static const std::array<Edges, TriadCensus::numberOfTriads> edges{{
            {},
            {{0, 1}},
            {{0, 1}, {1, 0}},
            {{1, 0}, {1, 2}},
            {{0, 1}, {2, 1}},
            {{0, 1}, {1, 2}},
            {{0, 1}, {1, 0}, {2, 1}},
            {{0, 1}, {1, 0}, {1, 2}},
            {{0, 1}, {2, 1}, {0, 2}},
            {{0, 1}, {1, 2}, {2, 0}},
            {{0, 1}, {1, 0}, {1, 2}, {2, 1}},
            {{1, 0}, {1, 2}, {0, 2}, {2, 0}},
            {{0, 1}, {2, 1}, {0, 2}, {2, 0}},
            {{0, 1}, {1, 2}, {0, 2}, {2, 0}},
            {{0, 1}, {1, 2}, {2, 1}, {0, 2}, {2, 0}},
            {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {0, 2}, {2, 0}},
        }};
        return edges;
    }

    static uint16_t maskOf(const Graph &G, node a, node b, node c) {
        const std::array<node, 3> nodes{a, b, c};
        uint16_t mask = 0;
        for (index i = 0; i < 3; ++i) {
            for (index j = 0; j < 3; ++j) {
                if (i != j && G.hasEdge(nodes[i], nodes[j])) {
                    mask |= uint16_t{1} << (i * 3 + j);
                }
            }
        }
        return mask;
    }

    static uint16_t maskOf(const Edges &edges, const std::array<node, 3> &labels) {
        uint16_t mask = 0;
        for (const auto &[u, v] : edges) {
            mask |= uint16_t{1} << (labels[u] * 3 + labels[v]);
        }
        return mask;
    }

    static std::vector<count> bruteForce(const Graph &G) {
        std::vector<uint8_t> triadOfMask(1 << 9, 0);
        for (index triad = 0; triad < TriadCensus::numberOfTriads; ++triad) {
            std::array<node, 3> labels{0, 1, 2};
            do {
                triadOfMask[maskOf(definitions()[triad], labels)] = static_cast<uint8_t>(triad);
            } while (std::next_permutation(labels.begin(), labels.end()));
        }

        std::vector<count> census(TriadCensus::numberOfTriads, 0);
        G.forNodes([&](node a) -> void {
            G.forNodes([&](node b) -> void {
                G.forNodes([&](node c) -> void {
                    if (a < b && b < c) {
                        ++census[triadOfMask[maskOf(G, a, b, c)]];
                    }
                });
            });
        });
        return census;
    }
};

TEST_F(TriadCensusGTest, testEachTriadUnderAllLabelings) {
    for (index triad = 0; triad < TriadCensus::numberOfTriads; ++triad) {
        std::array<node, 3> labels{0, 1, 2};
        do {
            Graph G(3, false, true);
            for (const auto &[u, v] : definitions()[triad]) {
                G.addEdge(labels[u], labels[v]);
            }

            TriadCensus census(G);
            census.run();

            std::vector<count> expected(TriadCensus::numberOfTriads, 0);
            expected[triad] = 1;
            EXPECT_EQ(census.getCounts(), expected) << "triad " << triad;
        } while (std::next_permutation(labels.begin(), labels.end()));
    }
}

TEST_F(TriadCensusGTest, testMatchesBruteForce) {
    for (const count seed : {1, 2, 3, 4, 5}) {
        Aux::Random::setSeed(seed, false);
        for (const double probability : {0.05, 0.2, 0.5, 0.9}) {
            constexpr count n = 25;
            Graph G(n, false, true);
            for (node u = 0; u < n; ++u) {
                for (node v = 0; v < n; ++v) {
                    if (u != v && Aux::Random::probability() < probability) {
                        G.addEdge(u, v);
                    }
                }
            }

            TriadCensus census(G);
            census.run();

            EXPECT_EQ(census.getCounts(), bruteForce(G));
        }
    }
}

TEST_F(TriadCensusGTest, testGetCount) {
    Graph G(4, false, true);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(0, 2);

    TriadCensus census(G);
    census.run();

    EXPECT_EQ(census.getCount(Triad::T030T), 1);
    EXPECT_EQ(census.getCount(Triad::T012), 3);
    EXPECT_EQ(census.getCount(Triad::T030C), 0);
}

TEST_F(TriadCensusGTest, testSelfLoopsAndParallelEdgesAreIgnored) {
    Graph G(3, false, true);
    G.addEdge(0, 1);
    G.addEdge(0, 1);
    G.addEdge(1, 1);
    G.addEdge(1, 2);

    TriadCensus census(G);
    census.run();

    EXPECT_EQ(census.getCount(Triad::T021C), 1);
}

TEST_F(TriadCensusGTest, testDeletedNodes) {
    Graph G(6, false, true);
    G.addEdge(1, 3);
    G.addEdge(3, 5);
    G.addEdge(0, 3);
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);

    TriadCensus census(G);
    census.run();

    EXPECT_EQ(census.getCount(Triad::T021C), 1);
    const std::vector<count> counts = census.getCounts();
    EXPECT_EQ(std::accumulate(counts.begin(), counts.end(), count{0}), 1);
}

TEST_F(TriadCensusGTest, testFewNodes) {
    for (const count n : {0, 1, 2}) {
        Graph G(n, false, true);
        if (n == 2) {
            G.addEdge(0, 1);
        }

        TriadCensus census(G);
        census.run();

        EXPECT_EQ(census.getCounts(), std::vector<count>(TriadCensus::numberOfTriads, 0));
    }
}

TEST_F(TriadCensusGTest, testNumberOfTriplesDoesNotOverflow) {
    Graph G(3000000, false, true);

    TriadCensus census(G);
    census.run();

    EXPECT_EQ(census.getCount(Triad::T003), 4499995500001000000ULL);
}

TEST_F(TriadCensusGTest, testRunTwice) {
    Graph G(3, false, true);
    G.addEdge(0, 1);

    TriadCensus census(G);
    census.run();
    census.run();

    EXPECT_EQ(census.getCount(Triad::T012), 1);
}

TEST_F(TriadCensusGTest, testInvalidInput) {
    EXPECT_THROW(TriadCensus{Graph(3)}, std::runtime_error);
}

TEST_F(TriadCensusGTest, testResultsBeforeRunThrow) {
    Graph G(3, false, true);
    TriadCensus census(G);
    EXPECT_THROW(census.getCount(Triad::T003), std::runtime_error);
    EXPECT_THROW(census.getCounts(), std::runtime_error);
}

} /* namespace NetworKit */

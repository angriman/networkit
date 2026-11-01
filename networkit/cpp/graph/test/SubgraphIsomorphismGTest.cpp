// networkit-format
#include <algorithm>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/SubgraphIsomorphism.hpp>

namespace NetworKit {

class SubgraphIsomorphismGTest : public testing::Test {
protected:
    static Graph triangle(bool directed = false) {
        Graph P(3, false, directed);
        P.addEdge(0, 1);
        P.addEdge(1, 2);
        P.addEdge(2, 0);
        return P;
    }

    static Graph path(count n, bool directed = false) {
        Graph P(n, false, directed);
        for (node u = 0; u + 1 < n; ++u) {
            P.addEdge(u, u + 1);
        }
        return P;
    }

    static Graph randomGraph(count n, double probability, bool directed) {
        Graph G(n, false, directed);
        for (node u = 0; u < n; ++u) {
            for (node v = 0; v < n; ++v) {
                if ((directed ? u != v : u < v) && Aux::Random::probability() < probability) {
                    G.addEdge(u, v);
                }
            }
        }
        return G;
    }

    static bool isMatch(const Graph &G, const Graph &P, const std::vector<node> &match,
                        bool induced) {
        bool valid = true;
        P.forNodePairs([&](node a, node b) -> void {
            if (match[a] == match[b]) {
                valid = false;
            }
            if (P.hasEdge(a, b) && !G.hasEdge(match[a], match[b])) {
                valid = false;
            }
            if (P.isDirected() && P.hasEdge(b, a) && !G.hasEdge(match[b], match[a])) {
                valid = false;
            }
            if (induced && G.hasEdge(match[a], match[b]) != P.hasEdge(a, b)) {
                valid = false;
            }
            if (induced && P.isDirected() && G.hasEdge(match[b], match[a]) != P.hasEdge(b, a)) {
                valid = false;
            }
        });
        return valid;
    }

    // All matches of a pattern with consecutive ids, trying all injective assignments.
    static std::vector<std::vector<node>> bruteForce(const Graph &G, const Graph &P, bool induced) {
        std::vector<std::vector<node>> result;
        std::vector<node> current(P.numberOfNodes());

        auto extend = [&](auto &&self, index position) -> void {
            if (position == current.size()) {
                if (isMatch(G, P, current, induced)) {
                    result.emplace_back(current.begin(), current.end());
                }
                return;
            }
            for (node u = 0; u < G.numberOfNodes(); ++u) {
                if (std::find(current.begin(), current.begin() + position, u)
                    == current.begin() + position) {
                    current[position] = u;
                    self(self, position + 1);
                }
            }
        };
        extend(extend, 0);
        return result;
    }

    static std::vector<std::vector<node>> sortedMatches(const SubgraphIsomorphism &algo) {
        auto matches = algo.getMatches();
        std::sort(matches.begin(), matches.end());
        return matches;
    }
};

TEST_F(SubgraphIsomorphismGTest, testTrianglesInCompleteGraph) {
    Graph G(4);
    G.forNodePairs([&](node u, node v) -> void { G.addEdge(u, v); });

    const Graph algoPattern = triangle();
    SubgraphIsomorphism algo(G, algoPattern);
    algo.run();

    // 4 triangles with 6 automorphisms each.
    EXPECT_EQ(algo.numberOfMatches(), 24);
}

TEST_F(SubgraphIsomorphismGTest, testInducedExcludesAdditionalEdges) {
    const Graph G = triangle();

    const Graph inducedPattern = path(3);
    SubgraphIsomorphism induced(G, inducedPattern, true);
    induced.run();
    EXPECT_EQ(induced.numberOfMatches(), 0);

    const Graph notInducedPattern = path(3);
    SubgraphIsomorphism notInduced(G, notInducedPattern, false);
    notInduced.run();
    EXPECT_EQ(notInduced.numberOfMatches(), 6);
}

TEST_F(SubgraphIsomorphismGTest, testDirectedRespectsOrientation) {
    Graph cycle = triangle(true);
    const Graph onCyclePattern = triangle(true);
    SubgraphIsomorphism onCycle(cycle, onCyclePattern);
    onCycle.run();
    EXPECT_EQ(onCycle.numberOfMatches(), 3);

    Graph acyclic(3, false, true);
    acyclic.addEdge(0, 1);
    acyclic.addEdge(1, 2);
    acyclic.addEdge(0, 2);
    const Graph onAcyclicPattern = triangle(true);
    SubgraphIsomorphism onAcyclic(acyclic, onAcyclicPattern);
    onAcyclic.run();
    EXPECT_EQ(onAcyclic.numberOfMatches(), 0);
}

TEST_F(SubgraphIsomorphismGTest, testDisconnectedPattern) {
    // Two isolated pattern nodes need any two different target nodes without an edge.
    Graph P(2);
    Graph G(3);
    G.addEdge(0, 1);

    SubgraphIsomorphism induced(G, P, true);
    induced.run();
    EXPECT_EQ(induced.numberOfMatches(), 4);

    SubgraphIsomorphism notInduced(G, P, false);
    notInduced.run();
    EXPECT_EQ(notInduced.numberOfMatches(), 6);
}

TEST_F(SubgraphIsomorphismGTest, testSelfLoops) {
    Graph P(2);
    P.addEdge(0, 0);
    P.addEdge(0, 1);

    Graph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 1);

    SubgraphIsomorphism withLoop(G, P);
    withLoop.run();
    EXPECT_EQ(sortedMatches(withLoop), std::vector<std::vector<node>>({{1, 0}, {1, 2}}));

    // An induced match of an edge without a self-loop cannot contain node 1, which has one.
    const Graph withoutLoopPattern = path(2);
    SubgraphIsomorphism withoutLoop(G, withoutLoopPattern);
    withoutLoop.run();
    EXPECT_EQ(withoutLoop.numberOfMatches(), 0);

    const Graph notInducedPattern = path(2);
    SubgraphIsomorphism notInduced(G, notInducedPattern, false);
    notInduced.run();
    EXPECT_EQ(notInduced.numberOfMatches(), 4);
}

TEST_F(SubgraphIsomorphismGTest, testParallelEdgesDoNotDuplicateMatches) {
    Graph G(2);
    G.addEdge(0, 1);
    G.addEdge(0, 1);

    Graph P(2);
    P.addEdge(0, 1);
    P.addEdge(0, 1);

    SubgraphIsomorphism algo(G, P);
    algo.run();
    EXPECT_EQ(algo.numberOfMatches(), 2);
}

TEST_F(SubgraphIsomorphismGTest, testDeletedNodes) {
    Graph G(5);
    G.addEdge(0, 2);
    G.addEdge(2, 4);
    G.addEdge(1, 2);
    G.removeNode(1);

    Graph P(4);
    P.addEdge(0, 3);
    P.addEdge(3, 2);
    P.removeNode(1);

    SubgraphIsomorphism algo(G, P);
    algo.run();

    EXPECT_EQ(sortedMatches(algo),
              std::vector<std::vector<node>>({{0, none, 4, 2}, {4, none, 0, 2}}));
}

TEST_F(SubgraphIsomorphismGTest, testPatternLargerThanGraph) {
    const Graph G = path(2);
    const Graph P = path(3);

    SubgraphIsomorphism algo(G, P);
    algo.run();
    EXPECT_EQ(algo.numberOfMatches(), 0);
}

TEST_F(SubgraphIsomorphismGTest, testMaxMatches) {
    Graph G(5);
    G.forNodePairs([&](node u, node v) -> void { G.addEdge(u, v); });
    const Graph P = triangle();

    SubgraphIsomorphism algo(G, P, true, 7);
    algo.run();

    EXPECT_EQ(algo.numberOfMatches(), 7);
    for (const auto &match : algo.getMatches()) {
        EXPECT_TRUE(isMatch(G, P, match, true));
    }
}

TEST_F(SubgraphIsomorphismGTest, testMatchesBruteForce) {
    for (const count seed : {1, 2, 3, 4, 5, 6, 7, 8}) {
        Aux::Random::setSeed(seed, false);
        for (const bool directed : {false, true}) {
            for (const bool induced : {false, true}) {
                const Graph G = randomGraph(7, 0.4, directed);
                const Graph P = randomGraph(4, 0.5, directed);

                SubgraphIsomorphism algo(G, P, induced);
                algo.run();

                auto expected = bruteForce(G, P, induced);
                std::sort(expected.begin(), expected.end());
                EXPECT_EQ(sortedMatches(algo), expected);
            }
        }
    }
}

TEST_F(SubgraphIsomorphismGTest, testRunTwice) {
    const Graph G = triangle();
    const Graph P = path(2);

    SubgraphIsomorphism algo(G, P);
    algo.run();
    const auto first = sortedMatches(algo);
    algo.run();

    EXPECT_EQ(first.size(), 6);
    EXPECT_EQ(sortedMatches(algo), first);
}

TEST_F(SubgraphIsomorphismGTest, testInvalidInput) {
    const Graph G = triangle();
    const Graph empty;
    const Graph directed = triangle(true);
    EXPECT_THROW(SubgraphIsomorphism(G, empty), std::invalid_argument);
    EXPECT_THROW(SubgraphIsomorphism(G, G, true, 0), std::invalid_argument);
    EXPECT_THROW(SubgraphIsomorphism(G, directed), std::invalid_argument);
}

TEST_F(SubgraphIsomorphismGTest, testResultsBeforeRunThrow) {
    const Graph G = triangle();
    SubgraphIsomorphism algo(G, G);
    EXPECT_THROW(algo.getMatches(), std::runtime_error);
    EXPECT_THROW(algo.numberOfMatches(), std::runtime_error);
}

} /* namespace NetworKit */

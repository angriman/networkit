/*
 * TopClosenessGTest.cpp
 *
 *  Created on: 17.11.2022
 *      Author: cls, Fabian Brandt-Tumescheit
 */

#include <iomanip>
#include <iostream>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/centrality/Closeness.hpp>
#include <networkit/centrality/HarmonicCloseness.hpp>
#include <networkit/centrality/TopCloseness.hpp>
#include <networkit/centrality/TopHarmonicCloseness.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class TopClosenessGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    bool useFirstHeu() const noexcept;
    bool useSecondHeu() const noexcept;
};

bool TopClosenessGTest::useFirstHeu() const noexcept {
    return GetParam().first;
}

bool TopClosenessGTest::useSecondHeu() const noexcept {
    return GetParam().second;
}

class TopHarmonicClosenessGTest : public testing::TestWithParam<bool> {
protected:
    bool useNBBound() const noexcept;
};

bool TopHarmonicClosenessGTest::useNBBound() const noexcept {
    return GetParam();
}

INSTANTIATE_TEST_SUITE_P(InstantiationName, TopClosenessGTest,
                         testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                         std::make_pair(false, true), std::make_pair(true, true)));

INSTANTIATE_TEST_SUITE_P(InstantiationName, TopHarmonicClosenessGTest,
                         testing::Values(false, true));

TEST_P(TopClosenessGTest, testTopCloseness) {
    constexpr count size = 400;
    constexpr count k = 10;

    for (bool isDirected : {false, true}) {
        Aux::Random::setSeed(42, false);
        const auto G1 = DorogovtsevMendesGenerator(size).generate();
        Graph G(G1, false, isDirected);

        Closeness cc(G1, true, ClosenessVariant::GENERALIZED);
        cc.run();
        auto exactScores = cc.scores();
        auto ranking = cc.ranking();
        TopCloseness topcc(G, k, useFirstHeu(), useSecondHeu());
        topcc.run();
        auto scores = topcc.topkScoresList();
        EXPECT_EQ(topcc.topkNodesList().size(), k);
        for (count i = 0; i < k; i++) {
            EXPECT_DOUBLE_EQ(ranking[i].second, scores[i]);
        }
    }
}

TEST_P(TopClosenessGTest, testTopClosenessWithNodeList) {
    METISGraphReader reader;
    Graph G = reader.read("input/lesmis.graph");
    constexpr count k = 10;
    const std::vector<node> nodeList{0, 1, 2, 3, 4, 5, 11, 26, 48, 64};

    // We expect complete TopCloseness to not contain nodes 0-5, while
    // restricted should. The first element for both complete and
    // restricted TopCC should be 11. {26, 48, 64} should also be present
    // in both results (but position might differ).
    TopCloseness topC(G, k, useFirstHeu(), useSecondHeu());
    topC.run();
    auto topCNodes = topC.topkNodesList();
    topC.restrictTopKComputationToNodes(nodeList);
    topC.run();
    auto topCRNodes = topC.topkNodesList();
    auto topCRScores = topC.topkScoresList();
    EXPECT_EQ(topCNodes[0], topCRNodes[0]);
    EXPECT_THAT(topCRNodes, testing::IsSupersetOf({0, 1, 2, 3, 4, 5, 26, 48, 64}));
    EXPECT_THAT(topCNodes, testing::IsSupersetOf({26, 48, 64}));
    EXPECT_THAT(topCNodes, testing::Not(testing::IsSupersetOf({0, 1, 2, 3, 4, 5})));
    EXPECT_TRUE(std::is_sorted(topCRScores.begin(), topCRScores.end(), std::greater<node>()));
}

TEST_P(TopHarmonicClosenessGTest, testTopHarmonicCloseness) {
    const double tol = 1e-6;
    constexpr count k = 3;

    /*
     * A clique of 6 nodes with a path of 6 nodes hanging off node 5: the clique nodes have the
     * highest closeness, so the search for the top-3 must prune the path. Undirected graphs also
     * get a separate edge and an isolated node; directed graphs also get the diamond
     * 12 -> {13, 14} -> 15 -> 16 feeding the cycle 17 -> 18 -> 19 -> 17, where bounds on the number
     * of reachable nodes are not tight.
     */
    for (bool isDirected : {false, true}) {
        for (bool isWeighted : {false, true}) {
            if (isWeighted && useNBBound()) {
                continue;
            }
            Graph G(isDirected ? 20 : 14, isWeighted, isDirected);
            const auto addEdge = [&](node u, node v) -> void {
                G.addEdge(u, v, 1 + (u + v) % 3);
                if (isDirected && u < 6 && v < 6) {
                    G.addEdge(v, u, 1 + (u + v) % 3);
                }
            };
            for (node u = 0; u < 6; ++u) {
                for (node v = u + 1; v < 6; ++v) {
                    addEdge(u, v);
                }
            }
            for (node u = 5; u < 11; ++u) {
                addEdge(u, u + 1);
            }
            if (isDirected) {
                for (const auto &[u, v] : std::to_array<std::pair<node, node>>(
                         {{12, 13}, {12, 14}, {13, 15}, {14, 15}, {15, 16}, {16, 17}, {17, 18},
                          {18, 19}, {19, 17}})) {
                    addEdge(u, v);
                }
            } else {
                addEdge(12, 13);
            }
            HarmonicCloseness cc(G, false);
            cc.run();
            const auto ranking = cc.ranking();
            TopHarmonicCloseness topcc(G, k, useNBBound());
            topcc.run();

            auto topkScores = topcc.topkScoresList();
            EXPECT_EQ(topcc.topkNodesList().size(), k);
            EXPECT_EQ(topkScores.size(), k);

            topkScores = topcc.topkScoresList(true);

            for (count i = 0; i < topkScores.size(); ++i) {
                EXPECT_NEAR(ranking[i].second, topkScores[i], tol);
            }
            for (count i = k; i < topkScores.size(); ++i) {
                EXPECT_NEAR(topkScores[i], topkScores[k - 1], tol);
            }
        }
    }
}

TEST_P(TopHarmonicClosenessGTest, testTopHarmonicClosenessWithNodeList) {
    METISGraphReader reader;
    Graph G = reader.read("input/lesmis.graph");
    constexpr count k = 10;
    const std::vector<node> nodeList{0, 1, 2, 3, 4, 5, 6, 11, 27, 48};

    // We expect complete TopHarmonicCloseness to not contain nodes 0-6,
    // while restricted should. {11, 27, 48} should be present in both
    // results (but position might differ).
    TopHarmonicCloseness topHC(G, k, useNBBound());
    topHC.run();
    auto topHCNodes = topHC.topkNodesList();
    topHC.restrictTopKComputationToNodes(nodeList);
    topHC.run();
    auto topHCRNodes = topHC.topkNodesList();
    auto topHCRScores = topHC.topkScoresList();
    EXPECT_THAT(topHCRNodes, testing::IsSupersetOf({0, 1, 2, 3, 4, 5, 6, 11, 27, 48}));
    EXPECT_THAT(topHCNodes, testing::IsSupersetOf({11, 27, 48}));
    EXPECT_THAT(topHCNodes, testing::Not(testing::IsSupersetOf({0, 1, 2, 3, 4, 5, 6})));
    EXPECT_TRUE(std::is_sorted(topHCRScores.begin(), topHCRScores.end(), std::greater<node>()));
}

} /* namespace NetworKit */

// networkit-format
/*
 * ClusteredRandomGraphGeneratorGTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class ClusteredRandomGraphGeneratorGTest : public testing::Test {};

TEST_F(ClusteredRandomGraphGeneratorGTest, testClusteredRandomGraphGenerator) {
    Aux::Random::setSeed(42, false);
    const count n = 100, c = 10;
    const double pin = 0.5, pout = 0.01;
    ClusteredRandomGraphGenerator gen(n, c, pin, pout);
    Graph G = gen.generate();
    Partition part = gen.getCommunities();
    count nCommunities = part.getSubsetIds().size();
    EXPECT_EQ(n, G.numberOfNodes());
    EXPECT_EQ(n, G.upperNodeIdBound());
    EXPECT_GE(nCommunities, 1);
    EXPECT_LE(nCommunities, c);
}

TEST_F(ClusteredRandomGraphGeneratorGTest,
       testClusteredRandomGraphGeneratorCompleteIsolatedCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 1, 0);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        { // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that each community is an isolated clique
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                const index uIdx = commPartition.subsetOf(u);
                EXPECT_EQ(G.degree(u), subsetSizes.at(uIdx) - 1);
                G.forNeighborsOf(u, [&](node v) { EXPECT_EQ(uIdx, commPartition.subsetOf(v)); });
            }
        }
    }
}

TEST_F(ClusteredRandomGraphGeneratorGTest, testClusteredRandomGraphGeneratorCompleteCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 1, 0.01);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        { // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that each community is a clique
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                const index uIdx = commPartition.subsetOf(u);
                EXPECT_GE(G.degree(u), subsetSizes.at(uIdx) - 1);
                count leftToVisit = community.size() - 1;
                G.forNeighborsOf(u, [&](node v) {
                    if (uIdx == commPartition.subsetOf(v))
                        --leftToVisit;
                });
                EXPECT_EQ(leftToVisit, 0);
            }
        }
    }
}

TEST_F(ClusteredRandomGraphGeneratorGTest,
       testClusteredRandomGraphGeneratorIndependentCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 0, 1);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        { // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that each community is an independent set, and the vertices in
        // each community are connected to all the vertices outside the community.
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                EXPECT_EQ(G.degree(u), n - community.size());
                G.forNeighborsOf(u, [&](node v) {
                    EXPECT_NE(commPartition.subsetOf(u), commPartition.subsetOf(v));
                });
            }
        }
    }
}

TEST_F(ClusteredRandomGraphGeneratorGTest, testClusteredRandomGraphGeneratorSparseCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 0.01, 1);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        { // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that the vertices in each community are connected to all the
        // vertices outside the community.
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                EXPECT_GE(G.degree(u), n - community.size());
                count leftToVisit = n - community.size();
                G.forNeighborsOf(u, [&](node v) {
                    if (commPartition.subsetOf(u) != commPartition.subsetOf(v))
                        --leftToVisit;
                });
                EXPECT_EQ(leftToVisit, 0);
            }
        }
    }
}

} /* namespace NetworKit */

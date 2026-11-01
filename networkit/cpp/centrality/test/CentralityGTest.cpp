// networkit-format
/*
 * CentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <typeinfo>
#include <unordered_set>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/centrality/ApproxElectricalCloseness.hpp>
#include <networkit/centrality/ApproxSpanningEdge.hpp>
#include <networkit/centrality/ComplexPaths.hpp>
#include <networkit/centrality/ForestCentrality.hpp>
#include <networkit/centrality/GedWalk.hpp>
#include <networkit/centrality/KPathCentrality.hpp>
#include <networkit/centrality/LocalPartitionCoverage.hpp>
#include <networkit/centrality/PermanenceCentrality.hpp>
#include <networkit/centrality/Sfigality.hpp>
#include <networkit/centrality/SpanningEdgeCentrality.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/structures/Partition.hpp>

#include "CentralityTestFixture.hpp"

namespace NetworKit {

class CentralityGTest : public DirectedWeightedGTest {};

INSTANTIATE_TEST_SUITE_P(InstantiationName, CentralityGTest, directedWeightedParams());

TEST_F(CentralityGTest, testLocalPartitionCoverage) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    Graph G(6);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    ClusteringGenerator clustGen;
    Partition P = clustGen.makeContinuousBalancedClustering(G, 3);

    LocalPartitionCoverage LPC(G, P);
    LPC.run();
    EXPECT_EQ(LPC.maximum(), 1);
    EXPECT_NEAR(LPC.centralization(), 0.29, 0.01);
}

TEST_F(CentralityGTest, testSfigality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    Graph G(6);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Sfigality sf(G);
    sf.run();
    auto res = sf.scores();
    vector<double> expectedRes = {1, 1, 0, 0.5, 0.5, 0};
    for (index i = 0; i < res.size(); ++i)
        EXPECT_EQ(res[i], expectedRes[i]);
    EXPECT_EQ(sf.maximum(), 1);
}

TEST_F(CentralityGTest, runKPathCentrality) {
    METISGraphReader reader;
    Graph G = reader.read("input/lesmis.graph");

    KPathCentrality centrality(G);
    centrality.run();
}

TEST_F(CentralityGTest, testSimplePermanence) {
    Graph G(15, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 0);
    G.addEdge(2, 3);
    node v = 4;
    node u = 5;
    G.addEdge(v, 0);
    G.addEdge(v, 1);
    G.addEdge(v, 2);
    G.addEdge(u, 3);
    G.addEdge(u, 2);
    G.addEdge(u, 0);
    G.addEdge(6, 7);
    G.addEdge(7, 8);
    G.addEdge(u, 6);
    G.addEdge(u, 7);
    G.addEdge(u, 8);
    G.addEdge(v, 6);
    G.addEdge(v, 7);
    G.addEdge(9, 10);
    G.addEdge(10, 11);
    G.addEdge(u, 9);
    G.addEdge(v, 10);
    G.addEdge(v, 11);
    G.addEdge(12, 13);
    G.addEdge(13, 14);
    G.addEdge(12, 14);
    G.addEdge(v, 12);
    G.addEdge(v, 14);

    Partition P(G.upperNodeIdBound());
    P.setUpperBound(4);
    P[0] = 0;
    P[1] = 0;
    P[2] = 0;
    P[3] = 0;
    P[v] = 0;
    P[u] = 0;
    P[6] = 1;
    P[7] = 1;
    P[8] = 1;
    P[9] = 2;
    P[10] = 2;
    P[11] = 2;
    P[12] = 3;
    P[13] = 3;
    P[14] = 3;

    ASSERT_EQ(9u, G.degree(v));
    ASSERT_EQ(7u, G.degree(u));

    PermanenceCentrality perm(G, P);
    perm.run();
    EXPECT_DOUBLE_EQ(2.0 / 3.0, perm.getIntraClustering(u));
    EXPECT_DOUBLE_EQ(1, perm.getIntraClustering(v));

    EXPECT_NEAR(-0.19048, perm.getPermanence(u), 0.0005);
    EXPECT_NEAR(0.167, perm.getPermanence(v), 0.0005);
}

TEST_F(CentralityGTest, testApproxSpanningEdge) {
    Aux::Random::setSeed(42, false);

    // Testing a graph that is too small yields approximation errors; however, if the input
    // graph has > 100 nodes, running the exact algorithm takes too long.
    // Therefore, in this test we use the approximation algorithm implemented in
    // SpanningEdgeCentrality as baseline, and check that the values are within a 2-epsilon
    // approximation.
    Graph G = ErdosRenyiGenerator(300, 0.1, false).generate();
    G.indexEdges();
    constexpr double eps = 0.1;

    ApproxSpanningEdge apx(G, eps);
    apx.run();
    SpanningEdgeCentrality se(G, eps);
    se.runParallelApproximation();
    auto apxScores = apx.scores();
    auto exactScores = se.scores();

    G.forEdges([&](node /*u*/, node /*v*/, edgeweight /*w*/, edgeid eid) {
        EXPECT_NEAR(apxScores[eid], exactScores[eid], 2 * eps);
    });
}

TEST_P(CentralityGTest, testGedWalk) {
    Aux::Random::setSeed(42, true);
    constexpr count k = 3;
    constexpr double epsilon = 0.01;
    auto g = ErdosRenyiGenerator(20, 0.1, isDirected()).generate();

    for (const auto bs : {GedWalk::BoundStrategy::GEOMETRIC, GedWalk::BoundStrategy::SPECTRAL}) {
        for (const auto gs : {GedWalk::GreedyStrategy::LAZY, GedWalk::GreedyStrategy::STOCHASTIC}) {
            GedWalk gedWalk(g, k, epsilon, -1.0, bs, gs);
            gedWalk.run();

            const auto apxScore = gedWalk.getApproximateScore();
            EXPECT_GE(apxScore, 0);
            const auto apxGroup = gedWalk.groupMaxGedWalk();
            EXPECT_EQ(std::unordered_set<node>(apxGroup.begin(), apxGroup.end()).size(), k);

            double maxScore = 0.;
            std::vector<node> group;
            std::vector<node> max_group;
            std::vector<bool> reference(g.numberOfNodes(), false);
            std::fill(reference.end() - k, reference.end(), true);

            do {
                group.clear();
                for (count i = 0; i < reference.size(); ++i) {
                    if (reference[i]) {
                        group.push_back(i);
                    }
                }

                const auto curScore = gedWalk.scoreOfGroup(group.begin(), group.end(), epsilon);
                if (curScore > maxScore) {
                    maxScore = curScore;
                    max_group = group;
                }
            } while (std::next_permutation(reference.begin(), reference.end()));

            EXPECT_GE(apxScore, (1. - 1. / std::exp(1)) * maxScore - epsilon);
        }
    }
}

TEST_F(CentralityGTest, testForestCentrality) {
    {
        // Directed graphs are not supported
        Graph G(10, true, true);
        G.addEdge(0, 1);
        EXPECT_THROW(ForestCentrality(G, 0), std::runtime_error);

        // Non-augmented graph
        G = GraphTools::toUndirected(G);
        G = GraphTools::toUnweighted(G);
        EXPECT_THROW(ForestCentrality(G, 0), std::runtime_error);

        // Augmented but non-compact graph
        node root = GraphTools::augmentGraph(G);
        G.removeNode(5);
        EXPECT_THROW(ForestCentrality(G, root), std::runtime_error);
    }

    Aux::Random::setSeed(42, true);
    auto G = HyperbolicGenerator(200).generate();
    const node root = GraphTools::augmentGraph(G);

    const double eps = 0.05;

    ForestCentrality fc(G, root, eps);
    fc.run();
    EXPECT_GT(fc.getNumberOfSamples(), 0);
    const auto apxDiag = fc.getDiagonal();

    ApproxElectricalCloseness ele(G);
    const auto diag = ele.computeExactDiagonal();
    G.forNodes([&](node u) { EXPECT_NEAR(apxDiag[u], diag[u], eps); });

    EXPECT_EQ(fc.scores().size(), G.upperNodeIdBound());
}

TEST_F(CentralityGTest, testComplexPathsSingleNode) {
    /* Graph:
     0    3
      \  / \
       2   5 - 6
      / \ /     \
     1   4 - - - 7 */
    Graph G(8);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(4, 7);
    G.addEdge(5, 6);
    G.addEdge(6, 7);
    Aux::Random::setSeed(42, false);

    count threshold = 2;
    ComplexPathAlgorithm CP(G, threshold, NetworKit::ComplexPathAlgorithm::singleNode, 0);
    CP.run();

    auto adopters = CP.getAdopters();
    // only the neighbour(2) of startnode(0) will adopt
    for (auto adoptee : adopters)
        EXPECT_TRUE(adoptee == 0 || adoptee == 2);

    auto resG = CP.getComplexGraph();
    EXPECT_EQ(resG.numberOfNodes(), G.numberOfNodes());
    EXPECT_EQ(resG.numberOfEdges(), 4);

    // the remaining subgraph only contains the edges that are connected to node 2, because the
    // adoption does not spread further
    resG.forEdges([&](node u, node v) { EXPECT_TRUE(u == 2 || v == 2); });
}

TEST_F(CentralityGTest, testComplexPathsAllNodes) {
    /* Graph:
     0    3
      \  / \
       2   5 - 6
      / \ /     \
     1   4 - - - 7 */
    Graph G(8);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(4, 7);
    G.addEdge(5, 6);
    G.addEdge(6, 7);
    Aux::Random::setSeed(42, false);

    count numberOfPaths = 2;
    ComplexPathAlgorithm CP(G, numberOfPaths, NetworKit::ComplexPathAlgorithm::allNodes);
    CP.run();

    auto normalizedProbability = CP.getPLci();
    // cannot ensure same values, therefore, checking size and type of result
    EXPECT_EQ(normalizedProbability.size(), G.numberOfNodes());
    for (index i = 0; i < normalizedProbability.size(); ++i) {
        EXPECT_EQ(typeid(normalizedProbability[i]), typeid(double));
    }
}

} /* namespace NetworKit */

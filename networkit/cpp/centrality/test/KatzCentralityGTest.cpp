// networkit-format
/*
 * KatzCentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <algorithm>
#include <utility>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/DynKatzCentrality.hpp>
#include <networkit/centrality/KatzCentrality.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>

namespace NetworKit {

class KatzCentralityGTest : public testing::Test {};

// TODO: replace by smaller graph
TEST_F(KatzCentralityGTest, testKatzCentralityDirected) {
    Aux::Log::setLogLevel("QUIET");
    const auto G = SNAPGraphReader{}.read("input/wiki-Vote.txt");
    KatzCentrality kc(G, 5e-4);
    kc.run();
    EXPECT_EQ(kc.ranking().front().first, 699);
}

TEST_F(KatzCentralityGTest, testKatzTopk) {
    Aux::Log::setLogLevel("QUIET");
    const auto G = METISGraphReader{}.read("input/caidaRouterLevel.graph");

    KatzCentrality exactAlgo(G, 0, 1.0);
    DynKatzCentrality topAlgo(G, 100);
    exactAlgo.run();
    topAlgo.run();

    // We cannot compare the ranking as the algorithms might return different
    // rankings for nodes that have equal/nearly equal scores. Instead,
    // epsilon-compare the exact scores of the i-th node and the expected i-th
    // node.
    auto exactRanking = exactAlgo.ranking();
    auto topRanking = topAlgo.ranking();
    for (count i = 0; i < std::min(G.numberOfNodes(), static_cast<count>(100)); i++)
        EXPECT_NEAR(exactAlgo.score(topRanking[i].first), exactRanking[i].second, 1e-6);
}

TEST_F(KatzCentralityGTest, testKatzDynamicAddition) {
    Graph G = METISGraphReader{}.read("input/caidaRouterLevel.graph");
    DynKatzCentrality kc(G, 100);
    DEBUG("start kc run");
    kc.run();
    DEBUG("finish kc");
    Aux::Random::setSeed(42, false);
    node u, v;
    do {
        u = GraphTools::randomNode(G);
        v = GraphTools::randomNode(G);
    } while (G.hasEdge(u, v));
    GraphEvent e(GraphEvent::EDGE_ADDITION, u, v, 1.0);
    kc.update(e);
    G.addEdge(u, v);
    DynKatzCentrality kc2(G, 100);
    kc2.run();
    const edgeweight tol = 1e-9;
    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        INFO("i = ", i);
        G.forNodes([&](node u) { EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u]); });
    }
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

TEST_F(KatzCentralityGTest, testKatzDynamicDeletion) {
    Graph G = METISGraphReader{}.read("input/caidaRouterLevel.graph");
    DynKatzCentrality kc(G, 100);
    DEBUG("start kc run");
    kc.run();
    DEBUG("finish kc");
    std::pair<node, node> p = GraphTools::randomEdge(G);
    node u = p.first;
    node v = p.second;
    INFO("Deleting edge ", u, ", ", v);
    GraphEvent e(GraphEvent::EDGE_REMOVAL, u, v, 1.0);
    G.removeEdge(u, v);
    kc.update(e);
    DynKatzCentrality kc2(G, 100);
    kc2.run();
    const edgeweight tol = 1e-9;
    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        INFO("i = ", i);
        G.forNodes([&](node u) {
            if (kc.nPaths[i][u] != kc2.nPaths[i][u]) {
                INFO("i = ", i, ", node ", u, ", dyn kc paths: ", kc.nPaths[i][u],
                     ", stat paths: ", kc2.nPaths[i][u]);
            }
            EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u]);
        });
    }
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

TEST_F(KatzCentralityGTest, testKatzDynamicBuilding) {
    Graph GIn = METISGraphReader{}.read("input/hep-th.graph");

    // Find a single max-degree node and add its edges to G.
    // (This guarantees that alpha is correct.)
    node maxNode = 0;
    GIn.forNodes([&](node u) {
        if (GIn.degree(u) > GIn.degree(maxNode))
            maxNode = u;
    });

    Graph G(GIn.upperNodeIdBound());

    GIn.forEdgesOf(maxNode, [&](node u, edgeweight) { G.addEdge(maxNode, u); });

    // Now run the algo. and add other some edges to check the correctness of
    // the dynamic part.
    DynKatzCentrality dynAlgo(G, 100);
    dynAlgo.run();

    count edgesProcessed = 0;
    GIn.forEdges([&](node u, node v) {
        if (u == maxNode || v == maxNode)
            return;
        if (edgesProcessed > 1000)
            return;
        GraphEvent e(GraphEvent::EDGE_ADDITION, u, v, 1.0);
        dynAlgo.update(e);
        G.addEdge(u, v);
        edgesProcessed++;
    });

    DynKatzCentrality topAlgo(G, 100);
    topAlgo.run();

    auto topRanking = topAlgo.ranking();
    auto dynRanking = dynAlgo.ranking();
    for (count i = 0; i < std::min(G.numberOfNodes(), count{100}); i++)
        EXPECT_FALSE(dynAlgo.areDistinguished(topRanking[i].first, dynRanking[i].first))
            << "Nodes " << topRanking[i].first << " and " << dynRanking[i].first
            << " should not be distinguished!";
}

TEST_F(KatzCentralityGTest, testKatzDirectedAddition) {
    // Same graph as in testCoreDecompositionDirected.
    count n = 16;
    Graph G(n, false, true);

    G.addEdge(2, 4);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 7);
    G.addEdge(6, 7);

    G.addEdge(6, 8);
    G.addEdge(6, 9);
    G.addEdge(6, 11);
    G.addEdge(7, 12);
    G.addEdge(8, 9);

    G.addEdge(8, 10);
    G.addEdge(8, 11);
    G.addEdge(8, 13);
    G.addEdge(9, 10);
    G.addEdge(9, 11);

    G.addEdge(9, 13);
    G.addEdge(10, 11);
    G.addEdge(10, 13);
    G.addEdge(10, 14);
    G.addEdge(11, 13);

    G.addEdge(11, 14);
    G.addEdge(12, 15);
    G.addEdge(13, 14);
    G.addEdge(14, 15);

    DynKatzCentrality kc(G, 5);
    kc.run();

    node u, v;
    Aux::Random::setSeed(42, false);
    do {
        u = GraphTools::randomNode(G);
        v = GraphTools::randomNode(G);
    } while (G.hasEdge(u, v));
    GraphEvent e(GraphEvent::EDGE_ADDITION, u, v, 1.0);
    kc.update(e);
    G.addEdge(u, v);

    DynKatzCentrality kc2(G, 5);
    kc2.run();

    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        G.forNodes([&](node u) { EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u]); });
    }
    const edgeweight tol = 1e-9;
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

TEST_F(KatzCentralityGTest, testKatzDirectedDeletion) {
    // Same graph as in testCoreDecompositionDirected.
    count n = 16;
    Graph G(n, false, true);

    G.addEdge(2, 4);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 7);
    G.addEdge(6, 7);

    G.addEdge(6, 8);
    G.addEdge(6, 9);
    G.addEdge(6, 11);
    G.addEdge(7, 12);
    G.addEdge(8, 9);

    G.addEdge(8, 10);
    G.addEdge(8, 11);
    G.addEdge(8, 13);
    G.addEdge(9, 10);
    G.addEdge(9, 11);

    G.addEdge(9, 13);
    G.addEdge(10, 11);
    G.addEdge(10, 13);
    G.addEdge(10, 14);
    G.addEdge(11, 13);

    G.addEdge(11, 14);
    G.addEdge(12, 15);
    G.addEdge(13, 14);
    G.addEdge(14, 15);

    DynKatzCentrality kc(G, 5);
    kc.run();

    Aux::Random::setSeed(42, false);
    std::pair<node, node> p = GraphTools::randomEdge(G);
    node u = p.first;
    node v = p.second;
    INFO("Removing ", u, " -> ", v);
    GraphEvent e(GraphEvent::EDGE_REMOVAL, u, v, 1.0);
    G.removeEdge(u, v);
    kc.update(e);

    DynKatzCentrality kc2(G, 5);
    kc2.run();

    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        G.forNodes([&](node u) {
            EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u])
                << i << "-length paths ending in " << u << " do not match!";
        });
    }
    const edgeweight tol = 1e-9;
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

} /* namespace NetworKit */

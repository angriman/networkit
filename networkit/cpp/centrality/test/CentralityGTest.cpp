/*
 * CentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <iomanip>
#include <iostream>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/centrality/ApproxCloseness.hpp>
#include <networkit/centrality/ApproxEffectiveResistance.hpp>
#include <networkit/centrality/ApproxGroupBetweenness.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/Closeness.hpp>
#include <networkit/centrality/CoreDecomposition.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/centrality/DynKatzCentrality.hpp>
#include <networkit/centrality/DynTopHarmonicCloseness.hpp>
#include <networkit/centrality/EigenvectorCentrality.hpp>
#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/centrality/GroupCloseness.hpp>
#include <networkit/centrality/GroupDegree.hpp>
#include <networkit/centrality/HarmonicCloseness.hpp>
#include <networkit/centrality/KPathCentrality.hpp>
#include <networkit/centrality/KadabraBetweenness.hpp>
#include <networkit/centrality/KatzCentrality.hpp>
#include <networkit/centrality/LaplacianCentrality.hpp>
#include <networkit/centrality/LocalClusteringCoefficient.hpp>
#include <networkit/centrality/PageRank.hpp>
#include <networkit/centrality/PermanenceCentrality.hpp>
#include <networkit/centrality/SpanningEdgeCentrality.hpp>
#include <networkit/centrality/TopCloseness.hpp>
#include <networkit/centrality/TopHarmonicCloseness.hpp>
#include <networkit/distance/CommuteTimeDistance.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/KONECTGraphReader.hpp>
#include <networkit/io/NetworkitBinaryReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/structures/Cover.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class CentralityGTest : public testing::Test {};

TEST_F(CentralityGTest, testApproxEffectiveResistance) {
    Aux::Random::setSeed(1, true);
    auto getBicTriangle = []() {
        auto G1 = ErdosRenyiGenerator(100, 0.1).generate();
        auto G2 = ErdosRenyiGenerator(100, 0.1).generate();
        auto G = ErdosRenyiGenerator(200, 0.1).generate();
        G.append(G1);
        G.append(G2);
        G.addEdge(199, 200);
        G.addEdge(199, 300);
        count newNodes = 3;
        while (newNodes--)
            G.addNode();
        G.addEdge(400, 401);
        G.addEdge(400, 402);
        G.addEdge(401, 402);
        G.addEdge(399, 400);

        return G;
    };

    INFO("Reading graph");
    auto G = NetworkitBinaryReader().read("/home/angriman/graphs/livemocha.nkb");
//    if (G.isDirected())
//        G = GraphTools::toUndirected(G);
    if (G.isWeighted())
        G = GraphTools::toUnweighted(G);
    G.removeMultiEdges();
    G.removeSelfLoops();
//    G = ConnectedComponents::extractLargestConnectedComponent(G, true);
//    G.removeSelfLoops();
//    G.removeMultiEdges();
//    ConnectedComponents cc(G);
//    cc.run();
//    if (cc.numberOfComponents() > 1)
//        throw std::runtime_error("Error: too many components");

    INFO("Done");
    static constexpr double eps = 0.1;
    static constexpr double tol = 1e-9;

    // Run approximation
    ApproxEffectiveResistance apx(G, eps);
    apx.init();
    apx.numberOfUSTs = 1000;
    apx.run();
    INFO(apx.timeToSample);
    INFO(apx.timeDFS);
    INFO(apx.timeToAggregate);
}


TEST_F(CentralityGTest, testBetweennessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    const double tol = 1e-3;
    EXPECT_NEAR(0.0, bc[0], tol);
    EXPECT_NEAR(0.0, bc[1], tol);
    EXPECT_NEAR(15.0, bc[2], tol);
    EXPECT_NEAR(3.0, bc[3], tol);
    EXPECT_NEAR(3.0, bc[4], tol);
    EXPECT_NEAR(1.0, bc[5], tol);
}

TEST_F(CentralityGTest, testBetweenness2Centrality) {
    /* Graph:
          0    3
          \  / \
          2    5
          /  \ /
          1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    const double tol = 1e-3;
    EXPECT_NEAR(0.0, bc[0], tol);
    EXPECT_NEAR(0.0, bc[1], tol);
    EXPECT_NEAR(15.0, bc[2], tol);
    EXPECT_NEAR(3.0, bc[3], tol);
    EXPECT_NEAR(3.0, bc[4], tol);
    EXPECT_NEAR(1.0, bc[5], tol);
}

TEST_F(CentralityGTest, runApproxBetweennessSmallGraph) {
    /* Graph:
     0    3
      \  / \
       2   5
      / \ /
     1   4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    double epsilon = 0.1; // error
    double delta = 0.1;   // confidence
    ApproxBetweenness centrality(G, epsilon, delta);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    ASSERT_LE(centrality.scores().size(), 1.0 / (epsilon * epsilon));

    DEBUG("scores: ", bc);
}

TEST_F(CentralityGTest, testBetweennessCentralityWeighted) {
    /* Graph:
     0    3   6
      \  / \ /
       2 -- 5
      /  \ / \
     1    4   7

     Edges in the upper row have weight 3,
     the edge in the middle row has weight 1.5,
     edges in the lower row have weight 2.
    */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    const double tol = 1e-3;
    EXPECT_NEAR(0.0, bc[0], tol);
    EXPECT_NEAR(0.0, bc[1], tol);
    EXPECT_NEAR(23.0, bc[2], tol);
    EXPECT_NEAR(0.0, bc[3], tol);
    EXPECT_NEAR(0.0, bc[4], tol);
    EXPECT_NEAR(23.0, bc[5], tol);
    EXPECT_NEAR(0.0, bc[6], tol);
    EXPECT_NEAR(0.0, bc[7], tol);
}

// TODO: replace by smaller graph
TEST_F(CentralityGTest, testKatzCentralityDirected) {
    SNAPGraphReader reader;
    Graph G = reader.read("input/wiki-Vote.txt");
    KatzCentrality kc(G);

    DEBUG("start kc run");
    kc.run();
    DEBUG("finish kc");
    std::vector<std::pair<node, double>> kc_ranking = kc.ranking();
    // std::vector<double> kc_scores = kc.scores();

    EXPECT_EQ(kc_ranking[0].first, 699);
}

TEST_F(CentralityGTest, testKatzTopk) {
    METISGraphReader reader;
    Graph G = reader.read("input/caidaRouterLevel.graph");

    // Compute max. degree to ensure that we use the same alpha for both algos.
    count maxDegree = 0;
    G.forNodes([&](node u) {
        if (G.degree(u) > maxDegree)
            maxDegree = G.degree(u);
    });

    KatzCentrality exactAlgo(G, 1.0 / (maxDegree + 1), 1.0);
    DynKatzCentrality topAlgo(G, 100);
    exactAlgo.run();
    topAlgo.run();

    // We cannot compare the ranking as the algorithms might return different
    // rankings for nodes that have equal/nearly equal scores. Instead,
    // epsilon-compare the exact scores of the i-th node and the expected i-th
    // node.
    auto exactRanking = exactAlgo.ranking();
    auto topRanking = topAlgo.ranking();
    for (count i = 0; i < std::min(G.numberOfNodes(), count{100}); i++)
        EXPECT_NEAR(exactAlgo.score(topRanking[i].first),
                    exactRanking[i].second, 1e-6);
}

TEST_F(CentralityGTest, testKatzDynamicAddition) {
    METISGraphReader reader;
    Graph G = reader.read("input/caidaRouterLevel.graph");
    DynKatzCentrality kc(G, 100);
    DEBUG("start kc run");
    kc.run();
    DEBUG("finish kc");
    node u, v;
    do {
        u = G.randomNode();
        v = G.randomNode();
    } while (G.hasEdge(u, v));
    GraphEvent e(GraphEvent::EDGE_ADDITION, u, v, 1.0);
    kc.update(e);
    G.addEdge(u, v);
    DynKatzCentrality kc2(G, 100);
    kc2.run();
    const edgeweight tol = 1e-9;
    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        INFO("i = ", i);
        G.forNodes([&](node u) {
            // if (kc.nPaths[i][u] != kc2.nPaths[i][u]) {
            //	 INFO("i = ", i, ", node ", u, ", dyn kc paths: ",
            // kc.nPaths[i][u], ",
            // stat paths: ", kc2.nPaths[i][u]);
            // }
            EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u]);
        });
    }
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

TEST_F(CentralityGTest, testKatzDynamicDeletion) {
    METISGraphReader reader;
    Graph G = reader.read("input/caidaRouterLevel.graph");
    DynKatzCentrality kc(G, 100);
    DEBUG("start kc run");
    kc.run();
    DEBUG("finish kc");
    std::pair<node, node> p = G.randomEdge();
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
                INFO("i = ", i, ", node ", u,
                     ", dyn kc paths: ", kc.nPaths[i][u],
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

TEST_F(CentralityGTest, testKatzDynamicBuilding) {
    METISGraphReader reader;
    Graph GIn = reader.read("input/hep-th.graph");

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
        EXPECT_FALSE(
            dynAlgo.areDistinguished(topRanking[i].first, dynRanking[i].first))
            << "Nodes " << topRanking[i].first << " and " << dynRanking[i].first
            << " should not be distinguished!";
}

TEST_F(CentralityGTest, testKatzDirectedAddition) {
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
        u = G.randomNode();
        v = G.randomNode();
    } while (G.hasEdge(u, v));
    GraphEvent e(GraphEvent::EDGE_ADDITION, u, v, 1.0);
    kc.update(e);
    G.addEdge(u, v);

    DynKatzCentrality kc2(G, 5);
    kc2.run();

    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        G.forNodes(
            [&](node u) { EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u]); });
    }
    const edgeweight tol = 1e-9;
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

TEST_F(CentralityGTest, testKatzDirectedDeletion) {
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
    std::pair<node, node> p = G.randomEdge();
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

// TODO: replace by smaller graph
TEST_F(CentralityGTest, testPageRankDirected) {
    SNAPGraphReader reader;
    Graph G = reader.read("input/wiki-Vote.txt");
    PageRank pr(G);

    DEBUG("start pr run");
    pr.run();
    DEBUG("finish pr");
    std::vector<std::pair<node, double>> pr_ranking = pr.ranking();

    const double tol = 1e-3;
    EXPECT_EQ(pr_ranking[0].first, 699);
    EXPECT_NEAR(pr_ranking[0].second, 0.00432, tol);
}

TEST_F(CentralityGTest, testEigenvectorCentrality) {
    /* Graph:
     0    3   6
      \  / \ /
       2 -- 5
      /  \ / \
     1    4   7

     Edges in the upper row have weight 3,
     the edge in the middle row has weight 1.5,
     edges in the lower row have weight 2.
    */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    EigenvectorCentrality centrality(G);
    centrality.run();
    std::vector<double> cen = centrality.scores();

    // computed with Matlab
    const double tol = 1e-4;
    EXPECT_NEAR(0.2254, fabs(cen[0]), tol);
    EXPECT_NEAR(0.1503, fabs(cen[1]), tol);
    EXPECT_NEAR(0.5290, fabs(cen[2]), tol);
    EXPECT_NEAR(0.4508, fabs(cen[3]), tol);
    EXPECT_NEAR(0.3006, fabs(cen[4]), tol);
    EXPECT_NEAR(0.5290, fabs(cen[5]), tol);
    EXPECT_NEAR(0.2254, fabs(cen[6]), tol);
    EXPECT_NEAR(0.1503, fabs(cen[7]), tol);
}

TEST_F(CentralityGTest, testPageRankCentrality) {
    /* Graph:
     0    3   6
      \  / \ /
       2 -- 5
      /  \ / \
     1    4   7

     Edges in the upper row have weight 3,
     the edge in the middle row has weight 1.5,
     edges in the lower row have weight 2.
    */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    double damp = 0.85;
    PageRank centrality(G, damp);
    centrality.run();
    std::vector<double> cen = centrality.scores();

    // compare to Matlab results
    const double tol = 1e-4;
    EXPECT_NEAR(0.0753, fabs(cen[0]), tol);
    EXPECT_NEAR(0.0565, fabs(cen[1]), tol);
    EXPECT_NEAR(0.2552, fabs(cen[2]), tol);
    EXPECT_NEAR(0.1319, fabs(cen[3]), tol);
    EXPECT_NEAR(0.0942, fabs(cen[4]), tol);
    EXPECT_NEAR(0.2552, fabs(cen[5]), tol);
    EXPECT_NEAR(0.0753, fabs(cen[6]), tol);
    EXPECT_NEAR(0.0565, fabs(cen[7]), tol);
}

TEST_F(CentralityGTest, benchSequentialBetweennessCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    Betweenness bc(G);
    bc.run();
    std::vector<std::pair<node, double>> ranking = bc.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchParallelBetweennessCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    Betweenness bc(G);
    bc.run();
    std::vector<std::pair<node, double>> ranking = bc.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchEigenvectorCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    EigenvectorCentrality cen(G);
    cen.run();
    std::vector<std::pair<node, double>> ranking = cen.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchPageRankCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    double damp = 0.85;
    PageRank cen(G, damp);
    cen.run();
    std::vector<std::pair<node, double>> ranking = cen.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, runEstimateBetweenness) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");

    EstimateBetweenness abc2(G, 100);
    abc2.run();

    DEBUG("approximated betweenness scores: ", abc2.scores());
}

// FIXME look out for tolerance limit in paper sample nodes
TEST_F(CentralityGTest, testApproxClosenessCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");

    ApproxCloseness acc(G, 453, 0, true);
    acc.run();

    std::vector<double> acc_scores = acc.scores();

    ASSERT_EQ(acc_scores.size(), 453);

    // compare sampled approx closeness values vs real closeness values
    // its a good compromise to not let the exact closness algorithm run
    ASSERT_NEAR(acc_scores[0], 0.416206, 0.000001);
    ASSERT_NEAR(acc_scores[10], 0.355906, 0.000001);
    ASSERT_NEAR(acc_scores[87], 0.420465, 0.000001);
    ASSERT_NEAR(acc_scores[121], 0.38865, 0.000001);
    ASSERT_NEAR(acc_scores[178], 0.4, 0.000001);
    ASSERT_NEAR(acc_scores[254], 0.397188, 0.000001);
    ASSERT_NEAR(acc_scores[307], 0.398238, 0.000001);
    ASSERT_NEAR(acc_scores[398], 0.37604, 0.000001);
    ASSERT_NEAR(acc_scores[406], 0.360734, 0.000001);
    ASSERT_NEAR(acc_scores[446], 0.396491, 0.000001);
}

TEST_F(CentralityGTest, testApproxClosenessCentralityOnToyGraph) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    ApproxCloseness acc(G, 6, 0.1, false);
    acc.run();
    std::vector<double> cc = acc.scores();

    double maximum = acc.maximum();

    const double tol = 0.2;
    EXPECT_NEAR(0.1, cc[0], tol);
    EXPECT_NEAR(0.1, cc[1], tol);
    EXPECT_NEAR(0.166667, cc[2], tol);
    EXPECT_NEAR(0.125, cc[3], tol);
    EXPECT_NEAR(0.125, cc[4], tol);
    EXPECT_NEAR(0.1, cc[5], tol);
    EXPECT_NEAR(0.2, maximum, tol);

    ApproxCloseness acc2(G, 4, 0.1, true);
    acc2.run();
    std::vector<double> cc2 = acc2.scores();

    double maximum2 = acc2.maximum();

    EXPECT_NEAR(0.5, cc2[0], tol);
    EXPECT_NEAR(0.5, cc2[1], tol);
    EXPECT_NEAR(0.833335, cc2[2], tol);
    EXPECT_NEAR(0.625, cc2[3], tol);
    EXPECT_NEAR(0.625, cc2[4], tol);
    EXPECT_NEAR(0.5, cc2[5], tol);
    EXPECT_NEAR(0.2, maximum2, tol);
}

TEST_F(CentralityGTest, testEdgeBetweennessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.indexEdges();

    Betweenness centrality(G, false, true);
    centrality.run();
    std::vector<double> bc = centrality.edgeScores();

    const double tol = 1e-3;
    EXPECT_NEAR(10.0, bc[0], tol);
    EXPECT_NEAR(10.0, bc[1], tol);
    EXPECT_NEAR(10.0, bc[2], tol);
    EXPECT_NEAR(10.0, bc[3], tol);
    EXPECT_NEAR(6.0, bc[4], tol);
    EXPECT_NEAR(6.0, bc[5], tol);
}

TEST_F(CentralityGTest, debugEdgeBetweennessCentrality) {
    auto path = "input/PGPgiantcompo.graph";
    METISGraphReader reader;
    Graph G = reader.read(path);
    G.indexEdges();

    Betweenness centrality(G, false, true);
    centrality.run();
    std::vector<double> bc = centrality.edgeScores();
}

TEST_F(CentralityGTest, testClosenessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Closeness centrality(G, false, ClosenessVariant::generalized);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    double maximum = centrality.maximum();

    const double tol = 1e-3;
    EXPECT_NEAR(0.1, bc[0], tol);
    EXPECT_NEAR(0.1, bc[1], tol);
    EXPECT_NEAR(0.166667, bc[2], tol);
    EXPECT_NEAR(0.125, bc[3], tol);
    EXPECT_NEAR(0.125, bc[4], tol);
    EXPECT_NEAR(0.1, bc[5], tol);
    EXPECT_NEAR(0.2, maximum, tol);
}

TEST_F(CentralityGTest, testClosenessCentralityWeighted) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);
    Graph Gw(n, true, false);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.forEdges([&](node u, node v) { Gw.addEdge(u, v, 1.); });

    Closeness centrality(G, false, ClosenessVariant::generalized);
    Closeness cw(Gw, false, ClosenessVariant::generalized);
    centrality.run();
    cw.run();

    auto bc = centrality.scores();
    auto bc_ = centrality.scores();

    const double tol = 1e-6;
    for (node u = 0; u < n; ++u) {
        EXPECT_NEAR(bc[u], bc_[u], tol);
    }
}

TEST_F(CentralityGTest, testClosenessCentralityDirected) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n, false, true);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Closeness centrality(G, true, ClosenessVariant::generalized);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    const double tol = 1e-6;
    EXPECT_NEAR(0.4, bc[0], tol);
    EXPECT_NEAR(0.4, bc[1], tol);
    EXPECT_NEAR(0.45, bc[2], tol);
    EXPECT_NEAR(0.2, bc[3], tol);
    EXPECT_NEAR(0.2, bc[4], tol);
    EXPECT_NEAR(0, bc[5], tol);
}

TEST_F(CentralityGTest, testHarmonicClosenessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    HarmonicCloseness centrality(G, false);
    centrality.run();
    std::vector<double> hc = centrality.scores();

    double maximum = centrality.maximum();

    const double tol = 1e-3;
    EXPECT_NEAR(2.833, hc[0], tol);
    EXPECT_NEAR(2.833, hc[1], tol);
    EXPECT_NEAR(4.5, hc[2], tol);
    EXPECT_NEAR(3.5, hc[3], tol);
    EXPECT_NEAR(3.5, hc[4], tol);
    EXPECT_NEAR(3.1667, hc[5], tol);
    EXPECT_NEAR(1, maximum, tol);
}

TEST_F(CentralityGTest, runKPathCentrality) {
    METISGraphReader reader;
    Graph G = reader.read("input/lesmis.graph");

    KPathCentrality centrality(G);
    centrality.run();
}

TEST_F(CentralityGTest, testCoreDecompositionSimple) {
    count n = 3;
    Graph G(n);
    G.addEdge(0, 1);

    CoreDecomposition coreDec(G);
    coreDec.run();
    std::vector<double> coreness = coreDec.scores();

    EXPECT_EQ(1u, coreness[0]) << "expected coreness";
    EXPECT_EQ(1u, coreness[1]) << "expected coreness";
    EXPECT_EQ(0u, coreness[2]) << "expected coreness";
}

TEST_F(CentralityGTest, testCoreDecomposition) {
    count n = 16;
    Graph G(n);

    // 	// create graph used in Baur et al. and network analysis lecture
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

    EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
    EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

    // compute core decomposition
    CoreDecomposition coreDec(G);
    coreDec.run();
    std::vector<double> coreness = coreDec.scores();
    // init cores
    // init shells
    Cover cores = coreDec.getCover();
    Partition shells = coreDec.getPartition();

    EXPECT_EQ(0u, coreness[0]) << "expected coreness";
    EXPECT_EQ(0u, coreness[1]) << "expected coreness";
    EXPECT_EQ(1u, coreness[2]) << "expected coreness";
    EXPECT_EQ(1u, coreness[3]) << "expected coreness";
    EXPECT_EQ(1u, coreness[4]) << "expected coreness";
    EXPECT_EQ(1u, coreness[5]) << "expected coreness";
    EXPECT_EQ(3u, coreness[6]) << "expected coreness";
    EXPECT_EQ(2u, coreness[7]) << "expected coreness";
    EXPECT_EQ(4u, coreness[8]) << "expected coreness";
    EXPECT_EQ(4u, coreness[9]) << "expected coreness";
    EXPECT_EQ(4u, coreness[10]) << "expected coreness";
    EXPECT_EQ(4u, coreness[11]) << "expected coreness";
    EXPECT_EQ(2u, coreness[12]) << "expected coreness";
    EXPECT_EQ(4u, coreness[13]) << "expected coreness";
    EXPECT_EQ(3u, coreness[14]) << "expected coreness";
    EXPECT_EQ(2u, coreness[15]) << "expected coreness";

    // for (index e = 0; e < n; e++) {
    // 	EXPECT_EQ(cores.contains(e), true);
    // 	EXPECT_EQ(shells.contains(e), true);
    // }
    // EXPECT_EQ(cores.get, coreness[15]) << "expected coreness";

    // test throw runtime error for self-loop in graph
    Graph H(2);
    H.addEdge(0, 1);
    H.addEdge(1, 1);
    EXPECT_ANY_THROW(CoreDecomposition CoreDec(H));
}

TEST_F(CentralityGTest, benchCoreDecompositionLocal) {
    METISGraphReader reader;
    std::vector<std::string> filenames = {"caidaRouterLevel", "wing",
                                          "astro-ph", "PGPgiantcompo"};

    for (auto f : filenames) {
        std::string filename("input/" + f + ".graph");
        DEBUG("about to read file ", filename);
        Graph G = reader.read(filename);
        G.removeSelfLoops();
        CoreDecomposition coreDec(G, false);
        Aux::Timer timer;
        timer.start();
        coreDec.run();
        timer.stop();
        INFO("Time for ParK of ", filename, ": ", timer.elapsedTag());

        CoreDecomposition coreDec2(G, true);
        timer.start();
        coreDec2.run();
        timer.stop();
        INFO("Time for bucket queue based k-core decomposition of ", filename,
             ": ", timer.elapsedTag());

        G.forNodes(
            [&](node u) { EXPECT_EQ(coreDec.score(u), coreDec2.score(u)); });
    }
}

TEST_F(CentralityGTest, testCoreDecompositionDirected) {
    count n = 16;
    Graph G(n, false, true);

    // 	// create graph used in Baur et al. and network analysis lecture
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

    EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
    EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

    // compute core decomposition
    CoreDecomposition coreDec(G);
    coreDec.run();
    std::vector<double> coreness = coreDec.scores();

    EXPECT_EQ(0u, coreness[0]) << "expected coreness";
    EXPECT_EQ(0u, coreness[1]) << "expected coreness";
    EXPECT_EQ(1u, coreness[2]) << "expected coreness";
    EXPECT_EQ(1u, coreness[3]) << "expected coreness";
    EXPECT_EQ(1u, coreness[4]) << "expected coreness";
    EXPECT_EQ(1u, coreness[5]) << "expected coreness";
    EXPECT_EQ(3u, coreness[6]) << "expected coreness";
    EXPECT_EQ(2u, coreness[7]) << "expected coreness";
    EXPECT_EQ(4u, coreness[8]) << "expected coreness";
    EXPECT_EQ(4u, coreness[9]) << "expected coreness";
    EXPECT_EQ(4u, coreness[10]) << "expected coreness";
    EXPECT_EQ(4u, coreness[11]) << "expected coreness";
    EXPECT_EQ(2u, coreness[12]) << "expected coreness";
    EXPECT_EQ(4u, coreness[13]) << "expected coreness";
    EXPECT_EQ(3u, coreness[14]) << "expected coreness";
    EXPECT_EQ(2u, coreness[15]) << "expected coreness";
}

TEST_F(CentralityGTest, testLocalClusteringCoefficientUndirected) {
    count n = 16;
    Graph G(n, false, false);

    // 	// create graph used in Baur et al. and network analysis lecture
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

    EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
    EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

    // compute core decomposition
    LocalClusteringCoefficient lcc(G);
    lcc.run();
    std::vector<double> lccScores = lcc.scores();
    std::vector<double> reference = {0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.5, 0.0,
                                     0.8, 0.8, 0.8, 0.6666666666666666,
                                     0.0, 0.8, 0.5, 0.0};

    EXPECT_EQ(reference, lccScores);

    LocalClusteringCoefficient lccTurbo(G, true);
    lccTurbo.run();
    EXPECT_EQ(reference, lccTurbo.scores());

    // test throw runtime error for self-loop in graph
    Graph H(2);
    H.addEdge(0, 1);
    H.addEdge(1, 1);
    EXPECT_ANY_THROW(LocalClusteringCoefficient lcc(H));
}

TEST_F(CentralityGTest, testLocalClusteringCoefficientUndirected2) {
    Graph G(6, false, false);
    G.addEdge(1, 0);
    G.addEdge(2, 0);
    G.addEdge(2, 1);
    G.addEdge(3, 2);
    G.addEdge(3, 0);
    G.addEdge(3, 1);
    G.addEdge(4, 2);
    G.addEdge(4, 0);
    G.addEdge(5, 3);
    G.addEdge(5, 4);
    G.addEdge(5, 1);
    LocalClusteringCoefficient lcc(G);
    lcc.run();
    std::vector<double> lccScores = lcc.scores();
    std::vector<double> reference = {0.6666666666666666, 0.6666666666666666,
                                     0.6666666666666666, 0.6666666666666666,
                                     0.3333333333333333, 0.3333333333333333};

    EXPECT_EQ(reference, lccScores);
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

TEST_F(CentralityGTest, testTopClosenessDirected) {
    count size = 400;
    count k = 10;
    Graph G1 = DorogovtsevMendesGenerator(size).generate();
    Graph G(G1.upperNodeIdBound(), false, true);
    G1.forEdges([&](node u, node v) {
        G.addEdge(u, v);
        G.addEdge(v, u);
    });
    Closeness cc(G1, true, ClosenessVariant::generalized);
    cc.run();
    TopCloseness topcc(G, k, true, true);
    topcc.run();
    const edgeweight tol = 1e-7;
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc.topkScoresList()[i], tol);
    }
    TopCloseness topcc2(G, k, true, false);
    topcc2.run();
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc2.topkScoresList()[i], tol);
    }
}

TEST_F(CentralityGTest, testTopClosenessUndirected) {
    count size = 400;
    count k = 10;
    Graph G1 = DorogovtsevMendesGenerator(size).generate();
    Graph G(G1.upperNodeIdBound(), false, false);
    G1.forEdges([&](node u, node v) {
        G.addEdge(u, v);
        G.addEdge(v, u);
    });
    Closeness cc(G1, true, ClosenessVariant::generalized);
    cc.run();
    TopCloseness topcc(G, k, true, true);
    topcc.run();
    const edgeweight tol = 1e-7;
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc.topkScoresList()[i], tol);
    }
    TopCloseness topcc2(G, k, true, false);
    topcc2.run();
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc2.topkScoresList()[i], tol);
    }
}

TEST_F(CentralityGTest, testTopHarmonicClosenessDirected) {
    count size = 400;
    count k = 10;
    Graph G1 = DorogovtsevMendesGenerator(size).generate();
    Graph G(G1.upperNodeIdBound(), false, true);
    G1.forEdges([&](node u, node v) {
        G.addEdge(u, v);
        G.addEdge(v, u);
    });
    HarmonicCloseness cc(G1, false);
    cc.run();
    TopHarmonicCloseness topcc(G, k);
    topcc.run();
    const edgeweight tol = 1e-7;
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc.topkScoresList()[i], tol);
    }
}

TEST_F(CentralityGTest, testTopHarmonicClosenessUndirected) {
    count size = 400;
    count k = 10;
    Graph G1 = DorogovtsevMendesGenerator(size).generate();
    Graph G(G1.upperNodeIdBound(), false, false);
    G1.forEdges([&](node u, node v) {
        G.addEdge(u, v);
        G.addEdge(v, u);
    });
    HarmonicCloseness cc(G1, false);
    cc.run();
    TopHarmonicCloseness topcc(G, k);
    topcc.run();
    const edgeweight tol = 1e-7;
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc.topkScoresList()[i], tol);
    }
    TopHarmonicCloseness topcc2(G, k);
    topcc2.run();
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc2.topkScoresList()[i], tol);
    }
}

TEST_F(CentralityGTest, testLaplacianCentrality) {
    // The graph structure and reference values for the scores are taken from
    // Qi et al., Laplacian centrality: A new centrality measure for weighted
    // networks.
    //
    // See
    // https://math.wvu.edu/~cqzhang/Publication-files/my-paper/INS-2012-Laplacian-W.pdf.
    Graph G(6, true);

    G.addEdge(0, 1, 4);
    G.addEdge(0, 2, 2);
    G.addEdge(1, 2, 1);
    G.addEdge(1, 3, 2);
    G.addEdge(1, 4, 2);
    G.addEdge(4, 5, 1);

    LaplacianCentrality lc(G);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_EQ(140, scores[0]);
    EXPECT_EQ(180, scores[1]);
    EXPECT_EQ(56, scores[2]);
    EXPECT_EQ(44, scores[3]);
    EXPECT_EQ(52, scores[4]);
    EXPECT_EQ(8, scores[5]);
}

TEST_F(CentralityGTest, testLaplacianCentralityNormalized) {
    Graph G(6, true);

    G.addEdge(0, 1, 4);
    G.addEdge(0, 2, 2);
    G.addEdge(1, 2, 1);
    G.addEdge(1, 3, 2);
    G.addEdge(1, 4, 2);
    G.addEdge(4, 5, 1);

    LaplacianCentrality lc(G, true);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_EQ(0.70, scores[0]);
    EXPECT_EQ(0.90, scores[1]);
    EXPECT_EQ(0.28, scores[2]);
    EXPECT_EQ(0.22, scores[3]);
    EXPECT_EQ(0.26, scores[4]);
    EXPECT_EQ(0.04, scores[5]);
}

TEST_F(CentralityGTest, testLaplacianCentralityUnweighted) {
    Graph G(6);

    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(1, 4);
    G.addEdge(4, 5);

    LaplacianCentrality lc(G);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_EQ(18, scores[0]);
    EXPECT_EQ(34, scores[1]);
    EXPECT_EQ(18, scores[2]);
    EXPECT_EQ(10, scores[3]);
    EXPECT_EQ(16, scores[4]);
    EXPECT_EQ(6, scores[5]);
}

TEST_F(CentralityGTest, testGroupDegreeUndirected) {
    Aux::Random::setSeed(42, false);
    count nodes = 12;
    Graph g = ErdosRenyiGenerator(nodes, 0.3, false).generate();
    count k = 5;

    GroupDegree gd(g, k, false);
    gd.run();
    count score = gd.getScore();
    GroupDegree gdIncludeGroup(g, k, true);
    gdIncludeGroup.run();
    count scorePlusGroup = gdIncludeGroup.getScore();

    std::vector<bool> reference(nodes, false);
    for (count i = nodes - k; i < nodes; ++i) {
        reference[i] = true;
    }

    auto computeGroupDegree = [&](std::vector<bool> curGroup, Graph g) {
        count result = 0;
        g.forNodes([&](node u) {
            if (!curGroup[u]) {
                bool neighborInGroup = false;
                g.forNeighborsOf(u, [&](node v) {
                    if (!neighborInGroup && curGroup[v]) {
                        neighborInGroup = true;
                        ++result;
                    }
                });
            }
        });
        return result;
    };

    count maxScore = 0;

    do {
        count curScore = computeGroupDegree(reference, g);
        if (curScore > maxScore) {
            maxScore = curScore;
        }
    } while (std::next_permutation(reference.begin(), reference.end()));

    EXPECT_TRUE(score > 0.5 * maxScore);
    EXPECT_TRUE(scorePlusGroup >
                (1.0 - 1.0 / std::exp(1.0) * (double)(maxScore + k)));
    EXPECT_EQ(score, gd.scoreOfGroup(gd.groupMaxDegree()));
    EXPECT_EQ(scorePlusGroup,
              gdIncludeGroup.scoreOfGroup(gdIncludeGroup.groupMaxDegree()));
}

TEST_F(CentralityGTest, testGroupDegreeDirected) {
    Aux::Random::setSeed(42, false);
    count nodes = 12;
    Graph g = ErdosRenyiGenerator(nodes, 0.3, true, false).generate();
    count k = 5;

    GroupDegree gd(g, k, false);
    gd.run();

    count scoreNoGroup = gd.getScore();
    GroupDegree gdIncludeGroup(g, k, true);
    gdIncludeGroup.run();
    count scorePlusGroup = gdIncludeGroup.getScore();

    std::vector<bool> reference(nodes, false);
    for (count i = nodes - k; i < nodes; ++i) {
        reference[i] = true;
    }

    auto computeGroupDegree = [&](std::vector<bool> curGroup, Graph g) {
        count result = 0;
        g.forNodes([&](node u) {
            if (!curGroup[u]) {
                bool neighborInGroup = false;
                g.forInNeighborsOf(u, [&](node v) {
                    if (!neighborInGroup && curGroup[v]) {
                        neighborInGroup = true;
                        ++result;
                    }
                });
            }
        });
        return result;
    };

    count maxScore = 0;

    do {
        count curScore = computeGroupDegree(reference, g);
        if (curScore > maxScore) {
            maxScore = curScore;
        }
    } while (std::next_permutation(reference.begin(), reference.end()));

    EXPECT_TRUE(scoreNoGroup > 0.5 * maxScore);
    EXPECT_TRUE(scorePlusGroup >
                (1.0 - 1.0 / std::exp(1.0)) * (double)(maxScore + k));
    EXPECT_EQ(scoreNoGroup, gd.scoreOfGroup(gd.groupMaxDegree()));
    EXPECT_EQ(scorePlusGroup,
              gdIncludeGroup.scoreOfGroup(gdIncludeGroup.groupMaxDegree()));
}

TEST_F(CentralityGTest, testGroupBetweennessScore) {
/** 0           5
 *  | \       / |
 *  1 - 3 - 4 - 6
 *  | /       \ |
 *  2           7
 */
    Graph G(8, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(0, 3);
    G.addEdge(1, 3);
    G.addEdge(2, 3);

    G.addEdge(3, 4);

    G.addEdge(4, 5);
    G.addEdge(4, 6);
    G.addEdge(4, 7);
    G.addEdge(5, 6);
    G.addEdge(6, 7);

    BFS bfs(G, 0, true, false);
    // Naively computes the group betweenness of a group of nodes S
    auto naiveGB = [&](const std::vector<node> &S) -> double {
        double score = 0;
        count n = G.upperNodeIdBound();
        std::vector<bool> inGroup(n);
        for (node u : S)
            inGroup[u] = true;
        for (node source = 0; source < n; ++source) {
            bfs.setSource(source);
            bfs.run();
            for (node target = 0; target < n; ++target) {
                if (target == source)
                    continue;

                auto paths = bfs.getPaths(target);
                if (paths.empty())
                    continue;
                double curScore = 0;
                for (auto &path : paths) {
                    for (node u : path) {
                        if (u != source && u != target && inGroup[u]) {
                            curScore += 1;
                            break;
                        }
                    }
                }

                score += curScore / static_cast<double>(paths.size());
            }
        }

        return score;
    };

    ApproxGroupBetweenness gb(G, 2, 0.1);
    count z = G.numberOfNodes();
    for (count k = 1; k <= z; ++k) {
        std::vector<bool> inGroup(z);
        for (index i = 0; i < k; ++i)
            inGroup[z - i - 1] = true;
        do {
            std::vector<node> group;
            group.reserve(k);
            for (node u = 0; u < z; ++u) {
                if (inGroup[u]) {
                    group.push_back(u);
                    if (group.size() == k)
                        break;
                }
            }
            double computedScore = gb.scoreOfGroup(group);
            double naiveScore = naiveGB(group);
            EXPECT_NEAR(computedScore, naiveScore, 1e-6);
        } while (std::next_permutation(inGroup.begin(), inGroup.end()));
    }
}

TEST_F(CentralityGTest, runTestApproxGroupBetweennessSmallGraph) {

    Aux::Random::setSeed(42, false);

    count n = 8;
    Graph g(n, false, false);

    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(0, 5);

    count k = 2;
    double eps = 0.1;
    ApproxGroupBetweenness gb(g, k, eps);
    gb.run();

    double maxScore = 0;
    std::vector<bool> inGroup(n);
    for (index i = 0; i < k; ++i)
        inGroup[n - i - 1] = true;
    do {
        std::vector<node> group;
        group.reserve(k);
        for (node u = 0; u < n; ++u) {
            if (inGroup[u]) {
                group.push_back(u);
                if (group.size() == k)
                    break;
            }
        }

        maxScore = std::max(maxScore, gb.scoreOfGroup(group));
    } while (std::next_permutation(inGroup.begin(), inGroup.end()));

    EXPECT_TRUE(gb.scoreOfGroup(gb.groupMaxBetweenness()) >= maxScore * eps);
}

TEST_F(CentralityGTest, testGroupCloseness) {
    Aux::Random::setSeed(42, false);

    Graph g(8, false, false);

    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(0, 5);

    count k = 3;

    GroupCloseness gc(g, k);
    gc.run();
    auto apx = gc.groupMaxCloseness();
    std::sort(apx.begin(), apx.end());
    std::vector<node> solution = {0, 2, 5};
    for (count i = 0; i < k; ++i) {
        EXPECT_EQ(apx[i], solution[i]);
    }

    EXPECT_NEAR(gc.scoreOfGroup(apx), 1.0, 1e-5);
}

/**
 * This test succeeds with the fixed random seed (42).
 * However, the Kadabra algorithm computes a correct epsilon-approximation of
 * the betweenness centrality score of all the nodes of the graph with high
 * probability. Thus, it is possible that, for a different random seed, this
 * test fails.
 */
TEST_F(CentralityGTest, testKadabraAbsolute) {
    Aux::Random::setSeed(42, true);
    const count n = 10;
    Graph g = ErdosRenyiGenerator(n, 0.1).generate();

    const double delta = 0.1;
    const double epsilon = 0.01;
    KadabraBetweenness kadabra(g, epsilon, delta);
    kadabra.run();
    auto scores = kadabra.topkScoresList();
    auto nodes = kadabra.topkNodesList();

    Betweenness betweenness(g, true);
    betweenness.run();
    count maxErrors = (count)std::ceil(delta * (double)n);

    count errors = 0;
    for (count i = 0; i < n; ++i) {
        if (std::abs(scores[i] - betweenness.score(nodes[i])) > delta) {
            ++errors;
        }
    }

    EXPECT_TRUE(errors <= maxErrors);
}

/**
 * This test succeeds with the fixed random seed (42).
 * However, the Kadabra algorithm finds the top-k nodes with
 * highest betweenness centrality with high probability. Thus, it is possible
 * that, for a different random seed, this test fails.
 */

TEST_F(CentralityGTest, testKadabraTopK) {
    Aux::Random::setSeed(42, true);
    const count n = 10;
    Graph g = ErdosRenyiGenerator(n, 0.1).generate();

    const double delta = 0.1;
    const double epsilon = 0.01;
    const count k = 3;
    KadabraBetweenness kadabra(g, epsilon, delta, k);
    kadabra.run();
    auto kadabraRanking = kadabra.ranking();

    Betweenness betweenness(g, true);
    betweenness.run();
    auto betwRanking = betweenness.ranking();
    bool correctRanking = true;
    for (count i = 0; i < k; ++i) {
        if (betwRanking[i].first != kadabraRanking[i].first) {
            correctRanking = false;
            int j = static_cast<int>(i) - 1;
            while (j >= 0 && betwRanking[j].second == betwRanking[i].second) {
                --j;
            }
            ++j;
            while (j < static_cast<int>(n) && betwRanking[j].second == betwRanking[i].second) {
                if (betwRanking[j].first == kadabraRanking[i].first) {
                    correctRanking = true;
                    break;
                }
                ++j;
            }
        }
    }
    EXPECT_TRUE(correctRanking);
}

TEST_F(CentralityGTest, testKadabraAbsoluteDeterministic) {
    const index seed = 42;
    Aux::Random::setSeed(seed, true);
    const count n = 300;
    Graph g = ErdosRenyiGenerator(n, 0.1).generate();

    const double delta = 0.1;
    const double epsilon = 0.01;
    Aux::Random::setSeed(seed, true);
    KadabraBetweenness kadabra(g, epsilon, delta, true);
    kadabra.run();
    auto scores = kadabra.topkScoresList();
    auto nodes = kadabra.topkNodesList();

    Aux::Random::setSeed(seed, true);
    KadabraBetweenness kadabra2(g, epsilon, delta, true);
    kadabra2.run();
    auto scores2 = kadabra2.topkScoresList();
    auto nodes2 = kadabra2.topkNodesList();

    for (count i = 0; i < n; ++i) {
        EXPECT_EQ(scores[i], scores2[i]);
    }
}

TEST_F(CentralityGTest, testDynTopHarmonicClosenessUndirected) {
    Graph G = DorogovtsevMendesGenerator(500).generate();

    count k = 10;

    DynTopHarmonicCloseness centrality(G, k, false);
    centrality.run();

    HarmonicCloseness reference(G, false);
    reference.run();

    auto scores = centrality.ranking();
    auto refScores = reference.ranking();

    for (count j = 0; j < k; ++j) {
        EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }

    count numInsertions = 1;

    std::vector<GraphEvent> deletions;
    std::vector<GraphEvent> insertions;

    for (count i = 0; i < numInsertions; i++) {

        node u = G.upperNodeIdBound();
        node v = G.upperNodeIdBound();

        do {
            u = G.randomNode();
            v = G.randomNode();
        } while (G.hasEdge(u, v));

        GraphEvent edgeAddition(GraphEvent::EDGE_ADDITION, u, v);
        insertions.insert(insertions.begin(), edgeAddition);

        GraphEvent edgeDeletion(GraphEvent::EDGE_REMOVAL, u, v);
        deletions.push_back(edgeDeletion);

        G.addEdge(u, v);
    }

    for (auto e : insertions) {
        G.removeEdge(e.u, e.v);
    }

    for (GraphEvent edgeAddition : insertions) {

        node u = edgeAddition.u;
        node v = edgeAddition.v;

        G.addEdge(u, v);

        centrality.update(edgeAddition);
        reference.run();

        scores = centrality.ranking();
        refScores = reference.ranking();

        for (count j = 0; j < k; ++j) {
            EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
        }
    }

    for (GraphEvent edgeDeletion : deletions) {
        node u = edgeDeletion.u;
        node v = edgeDeletion.v;

        G.removeEdge(u, v);

        centrality.update(edgeDeletion);
        reference.run();

        scores = centrality.ranking();
        refScores = reference.ranking();

        for (count j = 0; j < k; ++j) {
            EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
        }
    }

    for (GraphEvent edgeInsertion : insertions) {
        G.addEdge(edgeInsertion.u, edgeInsertion.v);
    }

    reference.run();
    centrality.updateBatch(insertions);

    scores = centrality.ranking();
    refScores = reference.ranking();

    for (count j = 0; j < k; ++j) {
        EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }
}

TEST_F(CentralityGTest, testDynTopHarmonicClosenessDirected) {
    Graph G = ErdosRenyiGenerator(300, 0.1, true).generate();

    count k = 10;

    DynTopHarmonicCloseness centrality(G, k, false);
    centrality.run();

    HarmonicCloseness reference(G, false);
    reference.run();

    auto scores = centrality.ranking();
    auto refScores = reference.ranking();
    for (count j = 0; j < k; ++j) {
        EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }

    count numInsertions = 1;

    std::vector<GraphEvent> deletions;
    std::vector<GraphEvent> insertions;

    for (count i = 0; i < numInsertions; i++) {

        node u = G.upperNodeIdBound();
        node v = G.upperNodeIdBound();

        do {
            u = G.randomNode();
            v = G.randomNode();
        } while (G.hasEdge(u, v));

        GraphEvent edgeAddition(GraphEvent::EDGE_ADDITION, u, v);
        insertions.insert(insertions.begin(), edgeAddition);

        GraphEvent edgeDeletion(GraphEvent::EDGE_REMOVAL, u, v);
        deletions.push_back(edgeDeletion);

        G.addEdge(u, v);
    }

    for (auto e : insertions) {
        G.removeEdge(e.u, e.v);
    }

    for (GraphEvent edgeAddition : insertions) {

        node u = edgeAddition.u;
        node v = edgeAddition.v;

        G.addEdge(u, v);

        centrality.update(edgeAddition);
        reference.run();

        scores = centrality.ranking();
        refScores = reference.ranking();

        for (count j = 0; j < k; ++j) {
            EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
        }
    }

    for (GraphEvent edgeDeletion : deletions) {

        node u = edgeDeletion.u;
        node v = edgeDeletion.v;

        G.removeEdge(u, v);

        centrality.update(edgeDeletion);
        reference.run();

        scores = centrality.ranking();
        refScores = reference.ranking();

        for (count j = 0; j < k; ++j) {
            EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
        }
    }

    for (GraphEvent edgeInsertion : insertions) {
        G.addEdge(edgeInsertion.u, edgeInsertion.v);
    }

    reference.run();
    centrality.updateBatch(insertions);

    scores = centrality.ranking();
    refScores = reference.ranking();

    for (count j = 0; j < k; ++j) {
        EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }
}

TEST_F(CentralityGTest, testSpielman) {
    SNAPGraphReader reader;
    auto G = reader.read("input/wiki-Vote.txt");
    G = ConnectedComponents::extractLargestConnectedComponent(G, true);
    SpanningEdgeCentrality sp(G);
    const auto diag = sp.computeDiagonal(0, 0.1, 1e-6);
    for (int i = 0; i < 10; ++i) {
        INFO(diag[i]);
    }
}

} /* namespace NetworKit */

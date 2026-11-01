// networkit-format
/*
 * GeneratorsGTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/dynamics/GraphUpdater.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/DynamicDorogovtsevMendesGenerator.hpp>
#include <networkit/generators/DynamicForestFireGenerator.hpp>
#include <networkit/generators/DynamicPathGenerator.hpp>
#include <networkit/generators/DynamicPubWebGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/MocnikGenerator.hpp>
#include <networkit/generators/MocnikGeneratorBasic.hpp>
#include <networkit/generators/PubWebGenerator.hpp>
#include <networkit/generators/RegularRingLatticeGenerator.hpp>
#include <networkit/generators/StochasticBlockmodel.hpp>
#include <networkit/generators/WattsStrogatzGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/viz/PostscriptWriter.hpp>

namespace NetworKit {

class GeneratorsGTest : public testing::Test {};

TEST_F(GeneratorsGTest, testStaticPubWebGenerator) {
    Aux::Random::setSeed(42, false);

    count n = 450;
    count numCluster = 9;
    count maxNumNeighbors = 36;
    float rad = 0.075;

    PubWebGenerator gen(n, numCluster, rad, maxNumNeighbors);
    Graph G = gen.generate();
    auto coordinates = gen.moveCoordinates();
    EXPECT_EQ(n, G.numberOfNodes()) << "number of generated nodes";

    // check degree
    G.forNodes([&](node v) { EXPECT_LE(G.degree(v), maxNumNeighbors) << "maximum degree"; });

    // 1-clustering
    ClusteringGenerator clusterGen;
    Partition oneClustering = clusterGen.makeOneClustering(G);
    EXPECT_EQ(G.numberOfNodes(), oneClustering.numberOfElements());

    // output to EPS file
    PostscriptWriter psWriter(true);
    psWriter.write(G, coordinates, oneClustering, "output/pubweb.eps");

    // clustering
    PLM clusterAlgo(G, false, 1.0, "none randomized");
    clusterAlgo.run();
    Partition clustering = clusterAlgo.getPartition();
    EXPECT_EQ(G.numberOfNodes(), clustering.numberOfElements());
    psWriter.write(G, coordinates, clustering, "output/pubweb-clustered-PLM.eps");

    Modularity mod;
    double modVal = mod.getQuality(clustering, G);
    EXPECT_GE(modVal, 0.2) << "modularity of clustering";
    DEBUG("Modularity of clustering: ", modVal);
    DEBUG("Total edge weight: ", G.totalEdgeWeight());
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testDynamicPubWebGenerator) {
    count nSteps = 5;
    count n = 200;
    count numCluster = 10;
    count maxNumNeighbors = 40;
    float rad = 0.08;

    DynamicPubWebGenerator dynGen(n, numCluster, rad, maxNumNeighbors, false);
    Graph G = dynGen.getGraph();
    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;

    for (index i = 1; i <= nSteps; ++i) {
        stream = dynGen.generate(1);
        DEBUG("updating graph");
        gu.update(stream);

        DEBUG("updated graph, new (n, m) = (", G.numberOfNodes(), ", ", G.numberOfEdges(), ")");
        edgeweight tew = G.totalEdgeWeight();
        DEBUG("1/2 graph volume: ", tew);
        EXPECT_GT(tew, 0);

// update coordinates
#ifndef NETWORKIT_RELEASE_LOGGING
        {
            PostscriptWriter psWriter(true);
            std::stringstream ss;
            ss << "output/pubweb-" << std::setw(4) << std::setfill('0') << i << ".eps";
            psWriter.write(G, dynGen.getCoordinates(), ss.str());
        }
#endif
    }
}

TEST_F(GeneratorsGTest, testDynamicPathGenerator) {
    count nSteps = 42;
    DynamicPathGenerator gen;
    auto stream = gen.generate(nSteps);
    EXPECT_EQ(stream.size(), nSteps * 3 + 1);
#ifndef NETWORKIT_RELEASE_LOGGING
    for (auto ev : stream) {
        TRACE(ev.toString());
    }
#endif
}

TEST_F(GeneratorsGTest, testErdosRenyiGenerator) {
    count n = 2000;
    double p = 1.5 * (std::log(n) / (double)n);

    ErdosRenyiGenerator generator(n, p);
    Graph G = generator.generate();
    EXPECT_EQ(n, G.numberOfNodes());
    EXPECT_FALSE(G.isEmpty());
    EXPECT_TRUE(G.checkConsistency());

    count nPairs = (n * (n - 1)) / 2;
    count nEdges = G.numberOfEdges();
    EXPECT_GE(nEdges, 0.75 * p * nPairs);
    EXPECT_LE(nEdges, 1.25 * p * nPairs);

    DEBUG("Number of edges with probability ", p, " (actual/expected): ", nEdges, " / ",
          (nPairs * p));
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testDynamicForestFireGenerator) {
    Graph G1(0);
    GraphUpdater gu1(G1);
    std::vector<GraphEvent> stream;
    DynamicForestFireGenerator ffg1(0.0, false);
    stream = ffg1.generate(10);
    gu1.update(stream);
    EXPECT_TRUE(G1.checkConsistency());
    EXPECT_EQ(10u, G1.numberOfNodes());
    G1.forNodes([&](node u) {
        count c = 0;
        G1.forNeighborsOf(u, [&](node v) {
            if (v < u) {
                c += 1;
            }
        });
        if (u == 0) {
            EXPECT_EQ(0u, c);
        } else {
            EXPECT_EQ(1u, c);
        }
    });

    Graph G2(0);
    GraphUpdater gu2(G2);
    DynamicForestFireGenerator ffg2(1.0, true, 1.0);
    stream = ffg2.generate(10);
    gu2.update(stream);
    EXPECT_TRUE(G2.checkConsistency());
    EXPECT_EQ(10u, G2.numberOfNodes());
    G2.forNodePairs([&](node u, node v) {
        if (v < u) {
            EXPECT_TRUE(G2.hasEdge(u, v));
        }
    });
    stream = ffg2.generate(10);
    gu2.update(stream);
    EXPECT_EQ(20u, G2.numberOfNodes());
}

TEST_F(GeneratorsGTest, testRegularRingLatticeGenerator) {
    int n0 = 10;
    int neighbors = 2;
    auto testRingLattice = [&](Graph G) {
        EXPECT_EQ(n0, (int)G.numberOfNodes());
        EXPECT_EQ(n0 * neighbors, (int)G.numberOfEdges());
        G.forNodePairs([&](node u, node v) {
            int diff = std::abs((int)u - (int)v);
            if (u != v && (diff <= neighbors || diff >= n0 - neighbors)) {
                EXPECT_TRUE(G.hasEdge(u, v));
            } else {
                EXPECT_FALSE(G.hasEdge(u, v));
            }
        });
    };

    RegularRingLatticeGenerator rrlg = RegularRingLatticeGenerator(n0, neighbors);
    testRingLattice(rrlg.generate());
}

TEST_F(GeneratorsGTest, testWattsStrogatzGenerator) {
    int n0 = 10;
    int neighbors = 2;
    auto testRingLattice = [&](Graph G) {
        G.forNodePairs([&](node u, node v) {
            int diff = std::abs((int)u - (int)v);
            if (u != v && (diff <= neighbors || diff >= n0 - neighbors)) {
                EXPECT_TRUE(G.hasEdge(u, v));
            } else {
                EXPECT_FALSE(G.hasEdge(u, v));
            }
        });
    };

    WattsStrogatzGenerator wsg1 = WattsStrogatzGenerator(n0, neighbors, 0.0);
    testRingLattice(wsg1.generate());

    WattsStrogatzGenerator wsg2 = WattsStrogatzGenerator(n0, neighbors, 0.3);
    Graph G = wsg2.generate();
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(n0, (int)G.numberOfNodes());
    EXPECT_EQ(n0 * neighbors, (int)G.numberOfEdges());
}

TEST_F(GeneratorsGTest, testWattsStrogatzGeneratorBigKs) {
    constexpr count nodes = 10;
    constexpr count neighbors = 4;
    constexpr double p = 0.5;
    for (int seed : {1, 2, 3, 4, 5}) {
        Aux::Random::setSeed(seed, false);
        const auto G = WattsStrogatzGenerator(nodes, neighbors, p).generate();
        EXPECT_TRUE(G.checkConsistency());
        EXPECT_EQ(nodes, G.numberOfNodes());
        EXPECT_EQ(nodes * neighbors, G.numberOfEdges());
    }
}

TEST_F(GeneratorsGTest, testDorogovtsevMendesGenerator) {
    int n0 = 20;
    DorogovtsevMendesGenerator dmg = DorogovtsevMendesGenerator(n0);
    Graph G = dmg.generate();

    EXPECT_EQ(n0, (int)G.numberOfNodes());
    EXPECT_EQ(2 * n0 - 3, (int)G.numberOfEdges());
    G.forNodes([&](node u) {
        count c = 0;
        G.forNeighborsOf(u, [&](node v) {
            if (v < u) {
                c += 1;
            }
        });
        if (u <= 2) {
            EXPECT_EQ(u, c);
        } else {
            EXPECT_EQ(2u, c);
        }
    });
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testDynamicDorogovtsevMendesGenerator) {
    count n0 = 20;
    DynamicDorogovtsevMendesGenerator ddmg = DynamicDorogovtsevMendesGenerator();
    Graph G(0);
    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;
    stream = ddmg.generate(n0 - 3);
    gu.update(stream);

    EXPECT_EQ(n0, G.numberOfNodes());
    EXPECT_EQ(2 * n0 - 3, G.numberOfEdges());
    G.forNodes([&](node u) {
        count c = 0;
        G.forNeighborsOf(u, [&](node v) {
            if (v < u) {
                c += 1;
            }
        });
        if (u <= 2) {
            EXPECT_EQ(u, c);
        } else {
            EXPECT_EQ(2u, c);
        }
    });
}

TEST_F(GeneratorsGTest, testStochasticBlockmodel) {
    count n = 10;
    count nBlocks = 2;
    std::vector<index> membership = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
    std::vector<std::vector<double>> affinity = {{1.0, 0.0}, {0.0, 1.0}};
    StochasticBlockmodel sbm(n, nBlocks, membership, affinity);
    Graph G = sbm.generate();

    EXPECT_EQ(n, G.numberOfNodes());
    EXPECT_EQ(20u, G.numberOfEdges());
}

TEST_F(GeneratorsGTest, testMocnikGenerator) {
    count dim = 3;
    count n = 10000;
    double k = 2.6;

    MocnikGenerator Mocnik(dim, n, k);
    Graph G(0);
    EXPECT_TRUE(G.isEmpty());
    G = Mocnik.generate();
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_NEAR(G.numberOfEdges() * 1. / G.numberOfNodes(), std::pow(k, dim), 20000);
}

TEST_F(GeneratorsGTest, testMocnikGeneratorBasic) {
    count dim = 3;
    count n = 5000;
    double k = 2.6;

    MocnikGeneratorBasic Mocnik(dim, n, k);
    Graph G(0);
    EXPECT_TRUE(G.isEmpty());
    G = Mocnik.generate();
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_NEAR(G.numberOfEdges() * 1. / G.numberOfNodes(), std::pow(k, dim), 10000);
}

} /* namespace NetworKit */

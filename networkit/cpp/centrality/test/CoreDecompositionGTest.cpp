// networkit-format
/*
 * CoreDecompositionGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/CoreDecomposition.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/structures/Cover.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class CoreDecompositionGTest : public testing::Test {};

TEST_F(CoreDecompositionGTest, testCoreDecompositionSimple) {
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

TEST_F(CoreDecompositionGTest, testCoreDecomposition) {
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

TEST_F(CoreDecompositionGTest, benchCoreDecompositionLocal) {
    METISGraphReader reader;
    std::vector<std::string> filenames = {"caidaRouterLevel", "wing", "astro-ph", "PGPgiantcompo"};

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

        CoreDecomposition coreDec2(G, false, true);
        timer.start();
        coreDec2.run();
        timer.stop();
        INFO("Time for bucket queue based k-core decomposition of ", filename, ": ",
             timer.elapsedTag());

        G.forNodes([&](node u) { EXPECT_EQ(coreDec.score(u), coreDec2.score(u)); });
    }
}

TEST_F(CoreDecompositionGTest, testCoreDecompositionDirected) {
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

} /* namespace NetworKit */

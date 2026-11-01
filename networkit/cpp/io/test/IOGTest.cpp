// networkit-format
/*
 * IOGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <tlx/unused.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/PLP.hpp>
#include <networkit/dynamics/GraphDifference.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/DGSReader.hpp>
#include <networkit/io/DotGraphWriter.hpp>
#include <networkit/io/EdgeListPartitionReader.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/EdgeListWriter.hpp>
#include <networkit/io/GraphIO.hpp>
#include <networkit/io/KONECTGraphReader.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/io/SNAPGraphWriter.hpp>
#include <networkit/io/ThrillGraphBinaryReader.hpp>
#include <networkit/io/ThrillGraphBinaryWriter.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class IOGTest : public testing::Test {};

TEST_F(IOGTest, testEdgeListWriter) {
    ErdosRenyiGenerator graphGen(100, 0.1, true);
    Graph G = graphGen.generate();

    std::string path = "output/edgelist2.txt";
    EdgeListWriter writer(' ', 1, true);
    EXPECT_NO_THROW(writer.write(G, path));

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    file.close();
    EXPECT_TRUE(exists) << "A file should have been created : " << path;

    EdgeListReader reader(' ', 1, "#", true, true);
    Graph G2 = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), G2.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), G2.numberOfEdges());
    EXPECT_EQ(G.isDirected(), G2.isDirected());
    EXPECT_EQ(G.isWeighted(), G2.isWeighted());

    // If not continuous, firstNode should be set to 0 automatically
    reader = EdgeListReader(' ', 1, "#", false, true);
    Graph G3 = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), G3.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), G3.numberOfEdges());
    EXPECT_EQ(G.isDirected(), G3.isDirected());
    EXPECT_EQ(G.isWeighted(), G3.isWeighted());
    std::vector<node> nodes(G.nodeRange().begin(), G.nodeRange().end());
    index i = 0;
    G3.forNodes([&](node u) { EXPECT_EQ(u, nodes[i++]); });
}

TEST_F(IOGTest, testGraphIOEdgeList) {
    ErdosRenyiGenerator graphGen(100, 0.1);
    Graph G = graphGen.generate();

    GraphIO graphio;
    std::string path = "output/edgelist.txt";
    graphio.writeEdgeList(G, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(IOGTest, testGraphIOAdjacencyList) {
    ErdosRenyiGenerator graphGen(100, 0.1);
    Graph G = graphGen.generate();
    GraphIO graphio;
    std::string path = "output/circular.adjlist";
    graphio.writeAdjacencyList(G, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(IOGTest, testGraphIOForIsolatedNodes) {
    Graph G(20);
    GraphIO graphio;
    std::string path = "output/isolated.adjlist";
    graphio.writeAdjacencyList(G, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(IOGTest, testDotGraphWriter) {
    ErdosRenyiGenerator graphGen(100, 0.1);
    Graph G = graphGen.generate();

    std::string path = "output/example.dot";

    DotGraphWriter writer;
    writer.write(G, path);

    // check if file exists
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "graph file should have been written to: " << path;
}

TEST_F(IOGTest, debugDGSReaderOnBigFile) {
    // read example graph
    DGSReader reader;
    Graph G;
    GraphEventProxy Gproxy(G);
    reader.read("/Users/forigem/KIT/NetworKit-CommunityDetection/input/AuthorsGraph.dgs", Gproxy);
}

TEST_F(IOGTest, debugDGSReader) {
    // read example graph
    DGSReader reader;
    Graph G;
    GraphEventProxy Gproxy(G);
    reader.read("input/example2.dgs", Gproxy);

    // get input parameters
    count nodeCount = G.numberOfNodes();
    DEBUG("Number of nodes ", nodeCount);
    EXPECT_EQ(3u, nodeCount);
    count edgeCount = G.numberOfEdges();
    DEBUG("Number of edges ", edgeCount);
    EXPECT_EQ(2u, edgeCount);

    G.forNodes([&](node n) {
        DEBUG("DEGREE OF NODE: ", G.degree(n), "\n");
        tlx::unused(n);
    });
}

TEST_F(IOGTest, testEdgeListReader) {
    EdgeListReader reader('\t', 1);

    std::string path = "input/network.dat";
    DEBUG("reading file: ", path);
    Graph G = reader.read(path);
    EXPECT_EQ(10u, G.numberOfNodes());
    EXPECT_EQ(10u, G.numberOfEdges());
    EXPECT_TRUE(G.hasEdge(0, 5));
    EXPECT_TRUE(G.hasEdge(2, 9));
    EXPECT_TRUE(G.hasEdge(1, 7));

    path = "input/example.edgelist";
    DEBUG("reading file: ", path);
    EdgeListReader reader2('\t', 1);
    Graph G2 = reader2.read(path);
    EXPECT_EQ(10u, G2.numberOfEdges());
    EXPECT_TRUE(G2.hasEdge(0, 4));

    path = "input/spaceseparated.edgelist";
    DEBUG("reading file: ", path);
    EdgeListReader reader3(' ', 1);
    Graph G3 = reader3.read(path);
    EXPECT_EQ(10u, G3.numberOfEdges());
    EXPECT_TRUE(G3.hasEdge(0, 4));

    path = "input/spaceseparated_weighted.edgelist";
    DEBUG("reading file: ", path);
    Graph G32 = reader3.read(path);
    EXPECT_TRUE(G32.isWeighted());
    EXPECT_EQ(2, G32.weight(0, 1));
    EXPECT_EQ(4, G32.weight(0, 2));
    EXPECT_EQ(3, G32.weight(1, 2));

    path = "input/comments.edgelist";
    DEBUG("reading file: ", path);
    EdgeListReader reader4('\t', 1);
    Graph G4 = reader4.read(path);
    EXPECT_EQ(10u, G4.numberOfEdges());
    EXPECT_TRUE(G4.hasEdge(0, 4));

    path = "input/alphabet.edgelist";
    DEBUG("reading file: ", path);
    EdgeListReader reader5('\t', 0, "#", false, false);
    Graph G5 = reader5.read(path);
    EXPECT_EQ(5u, G5.numberOfNodes());
    EXPECT_EQ(4u, G5.numberOfEdges());
    EXPECT_TRUE(G5.hasEdge(0, 1));
    EXPECT_TRUE(G5.hasEdge(0, 2));
    EXPECT_EQ(5, G5.weight(3, 4));
    EXPECT_EQ(1, G5.weight(2, 3));
}

TEST_F(IOGTest, debugReadingLFR) {
    std::string graphPath;
    std::string clustPath;

    std::cout << "[INPUT] LFR graph file path >" << std::endl;
    std::getline(std::cin, graphPath);

    std::cout << "[INPUT] clustering file path >" << std::endl;
    std::getline(std::cin, clustPath);

    EdgeListReader graphReader('\t', 1);
    EdgeListPartitionReader clusteringReader;

    Graph G = graphReader.read(graphPath);
    Partition truth = clusteringReader.read(clustPath);

    PLP PLP(G);
    PLP.run();
    Partition zeta = PLP.getPartition();

    Modularity mod;
    INFO("static clustering quality: ", mod.getQuality(zeta, G));
    INFO("static clustering number of clusters: ", zeta.numberOfSubsets());
    INFO("ground truth quality: ", mod.getQuality(truth, G));
    INFO("ground truth number of clusters: ", truth.numberOfSubsets());
}

TEST_F(IOGTest, debugReadingSNAP) {
    std::string graphPath;

    std::cout << "[INPUT] SNAP graph file path >" << std::endl;
    std::getline(std::cin, graphPath);

    EdgeListReader graphReader(' ', 1);

    Graph G = graphReader.read(graphPath);

    INFO("n = ", G.numberOfNodes());
    INFO("m = ", G.numberOfEdges());
}

TEST_F(IOGTest, testSNAPGraphReader) {
    SNAPGraphReader reader(true);
    Graph G = reader.read("input/wiki-Vote.txt");

    EXPECT_TRUE(G.isDirected());
    EXPECT_EQ(G.numberOfNodes(), 7115);
    EXPECT_EQ(G.numberOfEdges(), 103689);
}

TEST_F(IOGTest, testSNAPGraphWriter) {
    METISGraphReader reader;
    Graph G = reader.read("input/jazz.graph");

    std::string path = ""
                       "output/SNAPGraphWriter.gr";
    SNAPGraphWriter writer;
    writer.write(G, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "graph file should have been written to: " << path;
}

TEST_F(IOGTest, testThrillGraphBinaryWriterAndReader) {
    // This test graph has a large maximum degree as degrees smaller than 128
    // do not test the binary writer and reader properly.
    // So we simply use a star.
    count n = 257;
    Graph G(n, false, false);
    node center = 129;

    for (node u = 0; u < n; ++u) {
        if (u != center) {
            G.addEdge(u, center);
        }
    }

    std::string path = "output/test.thrillbin";

    ThrillGraphBinaryReader reader;
    ThrillGraphBinaryWriter writer;

    writer.write(G, path);
    Graph H = reader.read(path);

    GraphDifference diff(G, H);
    diff.run();

    EXPECT_EQ(diff.getEdits().size(), 0);
}

TEST_F(IOGTest, testKONECTGraphReader) {
    KONECTGraphReader reader;
    Graph G = reader.read("input/foodweb-baydry.konect");

    ASSERT_TRUE(G.isDirected());
    ASSERT_EQ(G.numberOfEdges(), 2137);
    ASSERT_EQ(G.numberOfNodes(), 128);
    ASSERT_EQ(G.weight(0, 1), 1.261404);
    ASSERT_EQ(G.weight(127, 48), 0.03050447);
}

} /* namespace NetworKit */

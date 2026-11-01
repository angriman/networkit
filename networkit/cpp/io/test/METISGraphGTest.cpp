// networkit-format
/*
 * METISGraphGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <fstream>
#include <iostream>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/METISGraphWriter.hpp>

namespace NetworKit {

class METISGraphGTest : public testing::Test {};

TEST_F(METISGraphGTest, testMETISGraphReader) {
    std::string path = "input/jazz.graph";

    METISGraphReader reader;
    Graph G = reader.read(path);
    count n = 198;
    count m = 2742;
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }

    path = "input/jazz2double.graph";
    G = reader.read(path);
    n = 5;
    m = 6;
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }

    // graph polblogs (has singletons)
    path = "input/polblogs.graph";
    G = reader.read(path);
    n = 1490;
    m = 16715;
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }

    // graph PGPgiantcompo
    path = "input/PGPgiantcompo.graph";
    G = reader.read(path);
    n = 10680;
    m = 24316;
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }
}

TEST_F(METISGraphGTest, testMETISGraphReaderWithTinyGraphs) {
    /* These graphs are from the METIS documentation and cover different settings
        of the fmt flag in the header of a METIS graph file */
    count n = 7;
    count m = 11;
    METISGraphReader reader;

    std::string path = "input/tiny_01.graph";
    Graph G = reader.read(path);
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }

    path = "input/tiny_02.graph";
    G = reader.read(path);
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }

    path = "input/tiny_03.graph";
    G = reader.read(path);
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }

    path = "input/tiny_04.graph";
    G = reader.read(path);
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }
}

TEST_F(METISGraphGTest, testMETISGraphWriter) {
    std::string path = "output/jazz1.graph";
    Graph G = Graph(3);
    G.addEdge(0, 2);
    G.addEdge(1, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 2);

    METISGraphWriter writer;
    writer.write(G, false, path);
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(METISGraphGTest, testMETISGraphWriterWithWeights) {
    std::string path = "output/jazz2.graph";
    Graph G = Graph(5);
    G.addEdge(0, 2);
    G.addEdge(0, 1);
    G.addEdge(0, 0);
    G.addEdge(1, 1);

    METISGraphWriter writer;
    writer.write(G, true, path);
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(METISGraphGTest, testMETISGraphReaderForNodeExistence2) {
    METISGraphReader reader;
    Graph G = reader.read("input/jazz.graph");
    EXPECT_TRUE(G.hasNode(0));
    EXPECT_EQ(198u, G.numberOfNodes());
    EXPECT_EQ(2742u, G.numberOfEdges());
}

TEST_F(METISGraphGTest, testMETISGraphReaderWithIsolatedNode) {
    METISGraphReader reader;
    Graph G = reader.read("input/example.graph");
    EXPECT_EQ(4u, G.numberOfNodes());
    EXPECT_EQ(2u, G.numberOfEdges());
    EXPECT_TRUE(G.hasNode(0));
    EXPECT_TRUE(G.hasNode(1));
    EXPECT_TRUE(G.hasNode(2));
    EXPECT_TRUE(G.hasNode(3));
    EXPECT_TRUE(G.hasEdge(0, 1));
    EXPECT_TRUE(G.hasEdge(0, 3));
}

TEST_F(METISGraphGTest, debugReadingMETISFile) {
    std::string graphPath;
    std::cout << "[INPUT] graph file path >" << std::endl;
    std::getline(std::cin, graphPath);

    METISGraphReader reader;
    Graph G = reader.read(graphPath);

    EXPECT_TRUE(true);
}

} /* namespace NetworKit */

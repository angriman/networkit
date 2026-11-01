// networkit-format
/*
 * GraphToolBinaryGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphToolBinaryReader.hpp>
#include <networkit/io/GraphToolBinaryWriter.hpp>

namespace NetworKit {

class GraphToolBinaryGTest : public testing::Test {};

TEST_F(GraphToolBinaryGTest, testGraphToolBinaryReader) {
    std::string path = "input/power.gt";
    GraphToolBinaryReader reader;
    Graph G = reader.read(path);
    EXPECT_EQ(4941u, G.numberOfNodes());
    EXPECT_EQ(6594u, G.numberOfEdges());
    EXPECT_FALSE(G.isDirected());
}

TEST_F(GraphToolBinaryGTest, testGraphToolBinaryWriter) {
    Graph G(10, false, false);
    G.addEdge(0, 1);
    G.addEdge(2, 1);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(5, 4);
    G.addEdge(5, 6);
    G.addEdge(7, 6);
    G.addEdge(8, 6);
    G.addEdge(7, 8);
    G.addEdge(9, 8);
    G.addEdge(9, 0);
    GraphToolBinaryReader reader;
    GraphToolBinaryWriter writer;
    std::string path = "output/test.gt";
    writer.write(G, path);
    Graph Gread = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), Gread.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), Gread.numberOfEdges());
    EXPECT_EQ(G.isDirected(), Gread.isDirected());
    EXPECT_EQ(G.isWeighted(), Gread.isWeighted());
}

TEST_F(GraphToolBinaryGTest, testGraphToolBinaryWriterDifferentGraphSizes) {
    GraphToolBinaryWriter writer;
    GraphToolBinaryReader reader;
    std::string path = "output/test.gt";
    // less than 8 nodes
    Graph G0(3, false, false);
    writer.write(G0, path);
    Graph G0read = reader.read(path);
    EXPECT_EQ(G0.numberOfNodes(), G0read.numberOfNodes());
    EXPECT_EQ(G0.numberOfEdges(), G0read.numberOfEdges());
    EXPECT_EQ(G0.isDirected(), G0read.isDirected());
    EXPECT_EQ(G0.isWeighted(), G0read.isWeighted());
    // less than 2^16 nodes
    Graph G1(30000, false, false);
    writer.write(G1, path);
    Graph G1read = reader.read(path);
    EXPECT_EQ(G1.numberOfNodes(), G1read.numberOfNodes());
    EXPECT_EQ(G1.numberOfEdges(), G1read.numberOfEdges());
    EXPECT_EQ(G1.isDirected(), G1read.isDirected());
    EXPECT_EQ(G1.isWeighted(), G1read.isWeighted());
    // less than 2^32 nodes
    Graph G2(100000, false, true);
    writer.write(G2, path);
    Graph G2read = reader.read(path);
    EXPECT_EQ(G2.numberOfNodes(), G2read.numberOfNodes());
    EXPECT_EQ(G2.numberOfEdges(), G2read.numberOfEdges());
    EXPECT_EQ(G2.isDirected(), G2read.isDirected());
    EXPECT_EQ(G2.isWeighted(), G2read.isWeighted());
}

TEST_F(GraphToolBinaryGTest, testGraphToolBinaryWriterWithDeletedNodes) {
    Graph G(10, false, false);
    G.removeNode(0);
    G.addEdge(2, 1);
    G.addEdge(2, 3);
    G.removeNode(4);
    G.addEdge(5, 6);
    G.addEdge(7, 6);
    G.addEdge(8, 6);
    G.addEdge(7, 8);
    G.removeNode(9);
    GraphToolBinaryReader reader;
    GraphToolBinaryWriter writer;
    std::string path = "output/test.gt";
    writer.write(G, path);
    Graph Gread = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), Gread.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), Gread.numberOfEdges());
    EXPECT_EQ(G.isDirected(), Gread.isDirected());
    EXPECT_EQ(G.isWeighted(), Gread.isWeighted());
}

TEST_F(GraphToolBinaryGTest, testGraphToolBinaryWriterDirected) {
    Graph G(10, false, true);
    G.addEdge(0, 1);
    G.addEdge(2, 1);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(5, 4);
    G.addEdge(5, 6);
    G.addEdge(7, 6);
    G.addEdge(8, 6);
    G.addEdge(7, 8);
    G.addEdge(9, 8);
    G.addEdge(9, 0);
    GraphToolBinaryReader reader;
    GraphToolBinaryWriter writer;
    std::string path = "output/test.gt";
    writer.write(G, path);
    Graph Gread = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), Gread.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), Gread.numberOfEdges());
    EXPECT_EQ(G.isDirected(), Gread.isDirected());
    EXPECT_EQ(G.isWeighted(), Gread.isWeighted());
}

TEST_F(GraphToolBinaryGTest, testGraphToolBinaryWriterWithDeletedNodesDirected) {
    Graph G(10, false, true);
    G.removeNode(0);
    G.addEdge(2, 1);
    G.addEdge(2, 3);
    G.removeNode(4);
    G.addEdge(5, 6);
    G.addEdge(7, 6);
    G.addEdge(8, 6);
    G.addEdge(7, 8);
    G.removeNode(9);
    GraphToolBinaryReader reader;
    GraphToolBinaryWriter writer;
    std::string path = "output/test.gt";
    writer.write(G, path);
    Graph Gread = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), Gread.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), Gread.numberOfEdges());
    EXPECT_EQ(G.isDirected(), Gread.isDirected());
    EXPECT_EQ(G.isWeighted(), Gread.isWeighted());
}

} /* namespace NetworKit */

// networkit-format
/*
 * GMLGraphGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <fstream>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GMLGraphReader.hpp>
#include <networkit/io/GMLGraphWriter.hpp>

namespace NetworKit {

class GMLGraphGTest : public testing::Test {};

TEST_F(GMLGraphGTest, testGMLGraphWriterUndirected) {
    std::string path = "output/jazz2_undirected.gml";
    Graph G = Graph(5);
    G.addEdge(0, 2);
    G.addEdge(0, 1);
    G.addEdge(0, 0);
    G.addEdge(1, 1);

    GMLGraphWriter writer;
    writer.write(G, path);
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(GMLGraphGTest, testGMLGraphWriterDirected) {
    std::string path = "output/jazz2_directed.gml";
    Graph G = Graph(5, false, true);
    G.addEdge(0, 2);
    G.addEdge(0, 1);
    G.addEdge(0, 0);
    G.addEdge(1, 1);

    GMLGraphWriter writer;
    writer.write(G, path);
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(GMLGraphGTest, testGMLGraphReaderUndirected) {
    std::string path = "input/jazz2_undirected.gml";
    GMLGraphReader reader;
    Graph G = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), 5u) << "number of nodes is not correct";
    EXPECT_TRUE(G.hasEdge(0, 2));
    EXPECT_TRUE(G.hasEdge(0, 1));
    EXPECT_TRUE(G.hasEdge(0, 0));
    EXPECT_TRUE(G.hasEdge(1, 1));
    EXPECT_FALSE(G.isDirected());
    EXPECT_TRUE(G.hasEdge(2, 0));
    EXPECT_TRUE(G.hasEdge(1, 0));
}

TEST_F(GMLGraphGTest, testGMLGraphReaderDirected) {
    std::string path = "input/jazz2_directed.gml";
    GMLGraphReader reader;
    Graph G = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), 5u) << "number of nodes is not correct";
    EXPECT_TRUE(G.hasEdge(0, 2));
    EXPECT_TRUE(G.hasEdge(0, 1));
    EXPECT_TRUE(G.hasEdge(0, 0));
    EXPECT_TRUE(G.hasEdge(1, 1));
    EXPECT_TRUE(G.isDirected());
    EXPECT_FALSE(G.hasEdge(2, 0));
    EXPECT_FALSE(G.hasEdge(1, 0));
}

} /* namespace NetworKit */

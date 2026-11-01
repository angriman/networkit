// networkit-format
/*
 * MatrixIOGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/MatrixTools.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/MTXGraphReader.hpp>
#include <networkit/io/MatrixMarketReader.hpp>
#include <networkit/io/RBMatrixReader.hpp>

namespace NetworKit {

class MatrixIOGTest : public testing::Test {};

TEST_F(MatrixIOGTest, testRBMatrixReader) {
    RBMatrixReader reader;
    auto csrRB = reader.read("input/tiny_05.rb");
    auto csrMM = MatrixMarketReader{}.read("input/tiny_05.matrixmarket");
    EXPECT_EQ(csrRB, csrMM);
}

TEST_F(MatrixIOGTest, testRBGraphReader) {
    RBGraphReader reader;
    auto graphRB = reader.read("input/tiny_05.rb");
    auto graphMM = MTXGraphReader{}.read("input/tiny_05.matrixmarket");

    EXPECT_EQ(graphRB.numberOfNodes(), graphMM.numberOfNodes());
    EXPECT_EQ(graphRB.numberOfEdges(), graphMM.numberOfEdges());
    EXPECT_EQ(graphRB.numberOfSelfLoops(), graphMM.numberOfSelfLoops());
    EXPECT_EQ(graphRB.isWeighted(), graphMM.isWeighted());
    EXPECT_EQ(graphRB.isDirected(), graphMM.isDirected());

    graphRB.forEdges([&](node u, node v, edgeweight w) {
        EXPECT_TRUE(graphMM.hasEdge(u, v));
        EXPECT_EQ(graphMM.weight(u, v), w);
    });

    graphMM.forEdges([&](node u, node v, edgeweight w) {
        EXPECT_TRUE(graphRB.hasEdge(u, v));
        EXPECT_EQ(graphRB.weight(u, v), w);
    });
}

TEST_F(MatrixIOGTest, testMatrixMarketReaderUnweightedUndirected) {
    CSRMatrix csr = MatrixMarketReader{}.read("input/chesapeake.mtx");
    EXPECT_EQ(csr.numberOfRows(), 39);
    EXPECT_EQ(csr.nnz(), 170 * 2); // double due to being undirected
    EXPECT_EQ(csr(7, 0), 1.0);
    EXPECT_EQ(csr(21, 1), 1.0);
    EXPECT_EQ(csr(5, 4), 1.0);
    EXPECT_EQ(csr(10, 6), 1.0);
    EXPECT_EQ(csr(35, 10), 1.0);
    EXPECT_EQ(csr(26, 13), 1.0);
    csr.forNonZeroElementsInRowOrder(
        [&](index i, index j, double) { EXPECT_EQ(csr(i, j), csr(j, i)); });

    // check that the MTXGraphReader is equivalent
    Graph G = MTXGraphReader{}.read("input/chesapeake.mtx");
    Graph G_from_csr = MatrixTools::matrixToGraph(csr);
    EXPECT_EQ(G.numberOfNodes(), G_from_csr.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), G_from_csr.numberOfEdges());
    G.forEdges([&](node u, node v, edgeweight w) {
        EXPECT_TRUE(G_from_csr.hasEdge(u, v));
        EXPECT_EQ(G_from_csr.weight(u, v), w);
    });
}

TEST_F(MatrixIOGTest, testMatrixMarketReaderUnweightedDirected) {
    CSRMatrix csr = MatrixMarketReader{}.read("input/GD01_b.mtx");
    EXPECT_EQ(csr.numberOfRows(), 18);
    EXPECT_EQ(csr.nnz(), 37);
    EXPECT_EQ(csr(7, 1), 1.0);
    EXPECT_EQ(csr(2, 3), 1.0);
    EXPECT_EQ(csr(15, 3), 1.0);
    EXPECT_EQ(csr(3, 9), 1.0);
    EXPECT_EQ(csr(11, 12), 1.0);
    EXPECT_EQ(csr(17, 13), 1.0);

    // check that the MTXGraphReader is equivalent
    Graph G = MTXGraphReader{}.read("input/GD01_b.mtx");
    Graph G_from_csr = MatrixTools::matrixToGraph(csr);
    EXPECT_EQ(G.numberOfNodes(), G_from_csr.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), G_from_csr.numberOfEdges());
    G.forEdges([&](node u, node v, edgeweight w) {
        EXPECT_TRUE(G_from_csr.hasEdge(u, v));
        EXPECT_EQ(G_from_csr.weight(u, v), w);
    });
}

TEST_F(MatrixIOGTest, testMatrixMarketReaderWeightedUndirected) {
    CSRMatrix csr = MatrixMarketReader{}.read("input/LFAT5.mtx");
    EXPECT_EQ(csr.numberOfRows(), 14);
    // double due to being undirected but we don't double count self loops
    EXPECT_EQ(csr.nnz(), 30 * 2 - 14);
    EXPECT_EQ(csr(5, 1), -6.2832e6);
    EXPECT_EQ(csr(3, 3), 15080.447999999997);
    EXPECT_EQ(csr(8, 3), 94.2528);
    EXPECT_EQ(csr(5, 5), 1.25664e7);
    EXPECT_EQ(csr(11, 8), -94.2528);
    EXPECT_EQ(csr(13, 11), 94.2528);
    csr.forNonZeroElementsInRowOrder(
        [&](index i, index j, double) { EXPECT_EQ(csr(i, j), csr(j, i)); });

    // check that the MTXGraphReader is equivalent
    Graph G = MTXGraphReader{}.read("input/LFAT5.mtx");
    Graph G_from_csr = MatrixTools::matrixToGraph(csr);
    EXPECT_EQ(G.numberOfNodes(), G_from_csr.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), G_from_csr.numberOfEdges());
    G.forEdges([&](node u, node v, edgeweight w) {
        EXPECT_TRUE(G_from_csr.hasEdge(u, v));
        EXPECT_EQ(G_from_csr.weight(u, v), w);
    });
}

TEST_F(MatrixIOGTest, testMatrixMarketReaderWeightedDirected) {
    CSRMatrix csr = MatrixMarketReader{}.read("input/Hamrle1.mtx");
    EXPECT_EQ(csr.numberOfRows(), 32);
    EXPECT_EQ(csr.nnz(), 98);
    EXPECT_EQ(csr(0, 0), 0.8499999999999999);
    EXPECT_EQ(csr(2, 2), 1.0);
    EXPECT_EQ(csr(3, 4), 0.2213493690543944);
    EXPECT_EQ(csr(5, 9), -0.2039265503510711);
    EXPECT_EQ(csr(9, 16), 0.2213493690543944);
    EXPECT_EQ(csr(11, 23), 0.04492168652740657);

    // check that the MTXGraphReader is equivalent
    Graph G = MTXGraphReader{}.read("input/Hamrle1.mtx");
    Graph G_from_csr = MatrixTools::matrixToGraph(csr);
    EXPECT_EQ(G.numberOfNodes(), G_from_csr.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), G_from_csr.numberOfEdges());
    G.forEdges([&](node u, node v, edgeweight w) {
        EXPECT_TRUE(G_from_csr.hasEdge(u, v));
        EXPECT_EQ(G_from_csr.weight(u, v), w);
    });
}

TEST_F(MatrixIOGTest, testMatrixMarketReaderIntegerWeights) {
    CSRMatrix csr = MatrixMarketReader{}.read("input/Ragusa16.mtx");
    EXPECT_EQ(csr.numberOfRows(), 24);
    EXPECT_EQ(csr.nnz(), 81);
    EXPECT_EQ(csr(0, 4), 1.0);
    EXPECT_EQ(csr(23, 4), 3.0);
    EXPECT_EQ(csr(10, 6), 1.0);
    EXPECT_EQ(csr(23, 9), 2.0);
    EXPECT_EQ(csr(4, 10), 4.0);
    EXPECT_EQ(csr(21, 21), 6.0);

    // check that the MTXGraphReader is equivalent
    Graph G = MTXGraphReader{}.read("input/Ragusa16.mtx");
    Graph G_from_csr = MatrixTools::matrixToGraph(csr);
    EXPECT_EQ(G.numberOfNodes(), G_from_csr.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), G_from_csr.numberOfEdges());
    G.forEdges([&](node u, node v, edgeweight w) {
        EXPECT_TRUE(G_from_csr.hasEdge(u, v));
        EXPECT_EQ(G_from_csr.weight(u, v), w);
    });
}

} /* namespace NetworKit */

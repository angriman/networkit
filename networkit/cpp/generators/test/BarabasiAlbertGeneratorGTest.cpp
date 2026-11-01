// networkit-format
/*
 * BarabasiAlbertGeneratorGTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <stdexcept>

#include <gtest/gtest.h>

#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/generators/DynamicBarabasiAlbertGenerator.hpp>
#include <networkit/generators/DynamicGraphSource.hpp>
#include <networkit/generators/PowerlawDegreeSequence.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphIO.hpp>
#include <networkit/io/METISGraphWriter.hpp>

namespace NetworKit {

class BarabasiAlbertGeneratorGTest : public testing::Test {};

TEST_F(BarabasiAlbertGeneratorGTest, testDynamicBarabasiAlbertGeneratorSingleStep) {
    count k = 2; // number of edges added per node
    DynamicGraphSource *gen = new DynamicBarabasiAlbertGenerator(k);
    GraphEventProxy *Gproxy = gen->newGraph();
    Graph *G = Gproxy->G;

    gen->initializeGraph();

    count nPre = G->numberOfNodes();
    count mPre = G->numberOfEdges();
    EXPECT_EQ(k, nPre) << "graph should have been initialized to k nodes";
    EXPECT_EQ(k - 1, mPre)
        << "graph should have been initialized to a path of k nodes which means k-1 edges";

    // perform single preferential attachment step
    gen->generate();

    count nPost = G->numberOfNodes();
    count mPost = G->numberOfEdges();
    EXPECT_EQ(nPre + 1, nPost) << "one more node should have been added";
    EXPECT_EQ(mPre + k, mPost) << "k edges should have been added";

    delete gen;
    delete Gproxy;
    delete G;
}

TEST_F(BarabasiAlbertGeneratorGTest, testDynamicBarabasiAlbertGenerator) {
    DynamicGraphSource *gen = new DynamicBarabasiAlbertGenerator(2);

    GraphEventProxy *Gproxy = gen->newGraph();
    Graph *G = Gproxy->G;

    gen->initializeGraph();

    EXPECT_EQ(2u, G->numberOfNodes()) << "initially the generator creates two connected nodes";
    EXPECT_EQ(1u, G->numberOfEdges()) << "initially the generator creates two connected nodes";

    count n = 100;

    gen->generateWhile([&]() { return (G->numberOfNodes() < n); });

    EXPECT_EQ(n, G->numberOfNodes());
    DEBUG("m = ", G->numberOfEdges());

    // resume generator

    gen->generateWhile([&]() { return (G->numberOfNodes() < 2 * n); });
    EXPECT_EQ(2 * n, G->numberOfNodes());

    delete gen;
    delete Gproxy;
    delete G;
}

TEST_F(BarabasiAlbertGeneratorGTest, viewDynamicBarabasiAlbertGenerator) {
    DynamicGraphSource *gen = new DynamicBarabasiAlbertGenerator(2);
    GraphEventProxy *Gproxy = gen->newGraph();
    Graph *G = Gproxy->G;
    gen->initializeGraph();
    count n = 42;
    gen->generateWhile([&]() { return (G->numberOfNodes() < n); });
    METISGraphWriter writer;
    writer.write(*G, "output/BATest.graph");

    delete gen;
    delete Gproxy;
    delete G;
}

TEST_F(BarabasiAlbertGeneratorGTest, testBarabasiAlbertGeneratorConstructor) {
    // k > nMax
    EXPECT_THROW(BarabasiAlbertGenerator generator(10, 9, 8, true), std::runtime_error);

    // n0 > nMax
    EXPECT_THROW(BarabasiAlbertGenerator generator(5, 9, 10, true), std::runtime_error);

    // n0 = initGraph.numberOfNodes() > nMax
    Graph initGraph(10);
    EXPECT_THROW(BarabasiAlbertGenerator generator(6, 9, initGraph, true), std::runtime_error);

    // initGraph, k > nMax
    initGraph = Graph(6);
    EXPECT_THROW(BarabasiAlbertGenerator generator(10, 9, initGraph, true), std::runtime_error);

    // initGraph does not have consecutive node ids
    initGraph.removeNode(0);
    EXPECT_THROW(BarabasiAlbertGenerator generator(3, 9, initGraph, true), std::runtime_error);
}

TEST_F(BarabasiAlbertGeneratorGTest, testBarabasiAlbertEmptyInit) {
    constexpr count k = 3;
    constexpr count nMax = 20;
    const Graph initGraph;

    BarabasiAlbertGenerator sequentialGenerator(k, nMax, initGraph, true);
    const Graph sequentialGraph = sequentialGenerator.generate();
    EXPECT_EQ(nMax, sequentialGraph.numberOfNodes());
    EXPECT_TRUE(sequentialGraph.checkConsistency());

    BarabasiAlbertGenerator parallelGenerator(k, nMax, initGraph, false);
    const Graph parallelGraph = parallelGenerator.generate();
    EXPECT_EQ(nMax, parallelGraph.numberOfNodes());
    EXPECT_TRUE(parallelGraph.checkConsistency());
}

TEST_F(BarabasiAlbertGeneratorGTest, testBarabasiAlbertParallelKOne) {
    constexpr count k = 1;
    constexpr count nMax = 8;
    const Graph initGraph;

    BarabasiAlbertGenerator generator(k, nMax, initGraph, false);
    const Graph G = generator.generate();

    EXPECT_EQ(nMax, G.numberOfNodes());
    EXPECT_TRUE(G.hasEdge(0, 1));
}

TEST_F(BarabasiAlbertGeneratorGTest, testBarabasiAlbertGeneratorBatagelj) {
    count k = 3;
    count nMax = 1000;
    count n0 = 3;

    BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0, true);
    Graph G = BarabasiAlbert.generate();

    EXPECT_EQ(nMax, G.numberOfNodes());
    EXPECT_LE(G.numberOfEdges(), nMax * k);
    EXPECT_TRUE(G.checkConsistency());

    Graph initGraph(4);
    initGraph.addEdge(0, 1);
    initGraph.addEdge(2, 1);
    initGraph.addEdge(2, 3);
    initGraph.addEdge(0, 3);
    BarabasiAlbert = BarabasiAlbertGenerator(k, nMax, initGraph, true);
    G = BarabasiAlbert.generate();

    EXPECT_EQ(nMax, G.numberOfNodes());
    EXPECT_LE(G.numberOfEdges(), nMax * k);
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(BarabasiAlbertGeneratorGTest, generatetBarabasiAlbertGeneratorGraph) {
    count k = 3;
    count nMax = 1000;
    count n0 = 3;

    BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);

    Graph G = BarabasiAlbert.generate();
    GraphIO io;
    io.writeAdjacencyList(G, "output/BarabasiGraph.txt");
}

TEST_F(BarabasiAlbertGeneratorGTest, testParallelBarabasiAlbertGeneratorDistribution) {
    count k = 3;
    count nMax = 1000;
    count n0 = 3;

    BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0, false);

    Graph G = BarabasiAlbert.generate();
    G.forNodes([&G, k, n0](node u) {
        if (u < n0)
            return;
        int totalEdges = 0;
        G.forEdgesOf(u, [&totalEdges, u](node v) {
            if (v < u)
                totalEdges++;
        });
        EXPECT_EQ(totalEdges, k);
    });
    PowerlawDegreeSequence sequence(G);
    auto gamma = sequence.getGamma();
    EXPECT_LT(gamma, -2);
    EXPECT_GT(gamma, -3);
}

} /* namespace NetworKit */

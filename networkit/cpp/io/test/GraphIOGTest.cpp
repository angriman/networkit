#include <cstdint>
#include <fstream>
#include <iterator>
#include <string>

#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/GraphIO.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using TestUtils::TemporaryDirectory;

template <class GraphT_>
class GraphIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;

protected:
    // Node ids from 10 on are printed as text also for 8-bit node types. Nodes 2 to 9 are
    // isolated.
    static GraphT smallGraph() {
        GraphT G(12);
        G.addEdge(NodeT{0}, NodeT{1});
        G.addEdge(NodeT{1}, NodeT{11});
        G.addEdge(NodeT{10}, NodeT{11});
        return G;
    }

    static std::string readFile(const std::string &path) {
        std::ifstream file(path);
        return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
    }
};

using GraphIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(GraphIOGTest, GraphIOTestTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(GraphIOGTest, testWriteEdgeList) {
    const TemporaryDirectory directory;
    const std::string path = directory.path("graph.edgelist");
    GraphIO{}.writeEdgeList(TestFixture::smallGraph(), path);
    EXPECT_EQ(TestFixture::readFile(path), "1 0\n"
                                           "11 1\n"
                                           "11 10\n");
}

TYPED_TEST(GraphIOGTest, testWriteAdjacencyList) {
    const TemporaryDirectory directory;
    const std::string path = directory.path("graph.adjlist");
    GraphIO{}.writeAdjacencyList(TestFixture::smallGraph(), path);
    EXPECT_EQ(TestFixture::readFile(path), "0 1\n"
                                           "1 0 11\n"
                                           "2\n3\n4\n5\n6\n7\n8\n9\n"
                                           "10 11\n"
                                           "11 1 10\n");
}

} // namespace
} // namespace NetworKit

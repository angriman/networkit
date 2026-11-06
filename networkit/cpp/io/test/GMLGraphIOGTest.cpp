#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/GMLGraphReader.hpp>
#include <networkit/io/GMLGraphWriter.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

template <class GraphT_>
class GMLGraphIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
};

using GMLGraphIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(GMLGraphIOGTest, GMLGraphIOTestTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(GMLGraphIOGTest, testRead) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    // The node ids of the file are labels, mapped to node ids in order of appearance.
    const std::string path = directory.writeFile("graph.gml", "graph [\n"
                                                              "  node [\n"
                                                              "    id b\n"
                                                              "  ]\n"
                                                              "  node [\n"
                                                              "    id a\n"
                                                              "  ]\n"
                                                              "  node [\n"
                                                              "    id c\n"
                                                              "  ]\n"
                                                              "  edge [\n"
                                                              "    source a\n"
                                                              "    target b\n"
                                                              "  ]\n"
                                                              "  edge [\n"
                                                              "    source a\n"
                                                              "    target c\n"
                                                              "  ]\n"
                                                              "]\n");
    const GraphT G = GenericGMLGraphReader<GraphT>{}.read(path);
    EXPECT_FALSE(G.isDirected());
    EXPECT_FALSE(G.isWeighted());
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{3});
    EXPECT_EQ(G.numberOfEdges(), count{2});
    EXPECT_TRUE(G.hasEdge(NodeT{1}, NodeT{0}));
    EXPECT_TRUE(G.hasEdge(NodeT{1}, NodeT{2}));
}

TYPED_TEST(GMLGraphIOGTest, testWriteAndReadBack) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;

    for (const bool directed : {false, true}) {
        // Node ids from 10 on are written as text also for 8-bit node types.
        GraphT G(12, false, directed);
        G.addEdge(NodeT{0}, NodeT{1});
        G.addEdge(NodeT{1}, NodeT{11});
        G.addEdge(NodeT{11}, NodeT{10});
        G.addEdge(NodeT{10}, NodeT{0});

        const std::string path = directory.path("roundtrip.gml");
        GenericGMLGraphWriter<GraphT>{}.write(G, path);
        const GraphT read = GenericGMLGraphReader<GraphT>{}.read(path);
        EXPECT_EQ(read.isDirected(), directed);
        EXPECT_TRUE(read.checkConsistency());
        EXPECT_EQ(read.numberOfNodes(), G.numberOfNodes());
        EXPECT_EQ(read.numberOfEdges(), G.numberOfEdges());
        G.forEdges([&](NodeT u, NodeT v) -> void { EXPECT_TRUE(read.hasEdge(u, v)); });
    }
}

TEST(GMLGraphIOGTest, testReadThrowsIfTooManyNodes) {
    using Reader = GenericGMLGraphReader<AdjListGraph<uint8_t, uint8_t>>;
    // The largest number of nodes that fits, 255, and then one more.
    constexpr uint32_t maxNodes = std::numeric_limits<uint8_t>::max();
    const auto gml = [](uint32_t n) -> std::string {
        std::string text = "graph [\n";
        for (uint32_t u = 0; u < n; ++u) {
            text += "  node [\n    id " + std::to_string(u) + "\n  ]\n";
        }
        return text + "]\n";
    };
    const TemporaryDirectory directory;
    EXPECT_EQ(Reader{}.read(directory.writeFile("fits.gml", gml(maxNodes))).numberOfNodes(),
              count{maxNodes});
    EXPECT_THAT(
        [&]() -> void { Reader{}.read(directory.writeFile("too_many.gml", gml(maxNodes + 1))); },
        ThrowsMessage<std::runtime_error>(HasSubstr("more nodes")));
}

} // namespace
} // namespace NetworKit

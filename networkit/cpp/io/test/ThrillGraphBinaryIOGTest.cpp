#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/ThrillGraphBinaryReader.hpp>
#include <networkit/io/ThrillGraphBinaryWriter.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

template <class GraphT_>
class ThrillGraphBinaryIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
};

using ThrillGraphBinaryIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(ThrillGraphBinaryIOGTest, ThrillGraphBinaryIOTestTypes,
                 /*Comma needed for variadic macro.*/);

TYPED_TEST(ThrillGraphBinaryIOGTest, testWriteAndReadBack) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;

    // A star, since degrees of at least 128 take more than one byte in the file. It has a self
    // loop and fewer than 255 nodes, so that it also fits 8-bit node types.
    constexpr count n = 200;
    const NodeT center{129};
    GraphT G(n);
    G.forNodes([&](NodeT u) -> void { G.addEdge(u, center); });

    const std::string path = directory.path("star.thrillbin");
    GenericThrillGraphBinaryWriter<GraphT>{}.write(G, path);

    for (const count nodesHint : {count{0}, n}) {
        const GraphT read = GenericThrillGraphBinaryReader<GraphT>(nodesHint).read(path);
        EXPECT_FALSE(read.isDirected());
        EXPECT_FALSE(read.isWeighted());
        EXPECT_TRUE(read.checkConsistency());
        EXPECT_EQ(read.numberOfNodes(), n);
        EXPECT_EQ(read.numberOfEdges(), n);
        EXPECT_EQ(read.numberOfSelfLoops(), count{1});
        EXPECT_EQ(read.degree(center), n);
        G.forEdges([&](NodeT u, NodeT v) -> void { EXPECT_TRUE(read.hasEdge(u, v)); });
    }
}

TYPED_TEST(ThrillGraphBinaryIOGTest, testReadMultipleFiles) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;

    // Node 0 is in the first file and node 1 in the second one.
    const std::string first =
        directory.writeFile("first.thrillbin", std::string("\x01\x01\x00\x00\x00", 5));
    const std::string second = directory.writeFile("second.thrillbin", std::string(1, '\0'));
    const GraphT G = GenericThrillGraphBinaryReader<GraphT>{}.read({first, second});
    EXPECT_EQ(G.numberOfNodes(), count{2});
    EXPECT_EQ(G.numberOfEdges(), count{1});
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{1}));
}

TEST(ThrillGraphBinaryIOGTest, testReadThrowsIfTooManyNodes) {
    using Reader = GenericThrillGraphBinaryReader<AdjListGraph<uint8_t, uint8_t>>;
    constexpr count maxNodes = std::numeric_limits<uint8_t>::max();
    const TemporaryDirectory directory;

    // One node per byte 0, which is an empty adjacency list.
    const std::string fits = directory.writeFile("fits.thrillbin", std::string(maxNodes, '\0'));
    EXPECT_EQ(Reader{}.read(fits).numberOfNodes(), maxNodes);
    const std::string tooMany =
        directory.writeFile("too_many.thrillbin", std::string(maxNodes + 1, '\0'));
    EXPECT_THAT([&]() -> void { Reader{}.read(tooMany); },
                ThrowsMessage<std::runtime_error>(HasSubstr("more nodes")));

    EXPECT_THAT([&]() -> void { Reader(maxNodes + 1).read(fits); },
                ThrowsMessage<std::runtime_error>(HasSubstr("does not fit")));
}

} // namespace
} // namespace NetworKit

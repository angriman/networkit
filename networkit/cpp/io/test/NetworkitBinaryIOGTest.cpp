#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/NetworkitBinaryReader.hpp>
#include <networkit/io/NetworkitBinaryWriter.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

template <class GraphT_>
class NetworkitBinaryIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

protected:
    // Returns a graph with a deleted node, a self loop, and weights representable by all tested
    // edge weight types.
    static GraphT smallGraph(bool weighted, bool directed) {
        GraphT G(13, weighted, directed);
        G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{2});
        G.addEdge(NodeT{1}, NodeT{12}, EdgeWeightT{5});
        G.addEdge(NodeT{12}, NodeT{11}, EdgeWeightT{100});
        G.addEdge(NodeT{11}, NodeT{0}, EdgeWeightT{7});
        G.addEdge(NodeT{11}, NodeT{11}, EdgeWeightT{3});
        G.removeNode(NodeT{5});
        return G;
    }

    // Checks that read has the same nodes, edges, weights, and edge ids as G.
    static void expectSameGraph(const GraphT &G, const GraphT &read) {
        EXPECT_EQ(read.isDirected(), G.isDirected());
        EXPECT_EQ(read.isWeighted(), G.isWeighted());
        EXPECT_EQ(read.hasEdgeIds(), G.hasEdgeIds());
        EXPECT_TRUE(read.checkConsistency());
        EXPECT_EQ(read.numberOfNodes(), G.numberOfNodes());
        EXPECT_EQ(read.upperNodeIdBound(), G.upperNodeIdBound());
        EXPECT_EQ(read.numberOfEdges(), G.numberOfEdges());
        EXPECT_EQ(read.numberOfSelfLoops(), G.numberOfSelfLoops());
        EXPECT_FALSE(read.hasNode(NodeT{5}));
        G.forEdges([&](NodeT u, NodeT v, EdgeWeightT w, edgeid id) -> void {
            EXPECT_TRUE(read.hasEdge(u, v));
            EXPECT_EQ(read.weight(u, v), w);
            if (G.hasEdgeIds()) {
                EXPECT_EQ(read.edgeId(u, v), id);
            }
        });
        if (G.hasEdgeIds()) {
            EXPECT_EQ(read.upperEdgeIdBound(), G.upperEdgeIdBound());
        }
    }
};

using NetworkitBinaryIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(NetworkitBinaryIOGTest, NetworkitBinaryIOTestTypes,
                 /*Comma needed for variadic macro.*/);

TYPED_TEST(NetworkitBinaryIOGTest, testWriteAndReadBack) {
    using GraphT = typename TestFixture::GraphT;
    const TemporaryDirectory directory;
    const std::string path = directory.path("graph.nkbg");

    for (const bool weighted : {false, true}) {
        for (const bool directed : {false, true}) {
            for (const bool indexed : {false, true}) {
                GraphT G = TestFixture::smallGraph(weighted, directed);
                if (indexed) {
                    G.indexEdges();
                }
                GenericNetworkitBinaryWriter<GraphT>{}.write(G, path);
                TestFixture::expectSameGraph(G, GenericNetworkitBinaryReader<GraphT>{}.read(path));

                const std::vector<uint8_t> buffer =
                    GenericNetworkitBinaryWriter<GraphT>{}.writeToBuffer(G);
                TestFixture::expectSameGraph(
                    G, GenericNetworkitBinaryReader<GraphT>{}.readFromBuffer(buffer));
            }
        }
    }
}

TYPED_TEST(NetworkitBinaryIOGTest, testWriteAndReadBackWithWeightFormats) {
    using GraphT = typename TestFixture::GraphT;
    for (const auto format :
         {NetworkitBinaryWeights::UNSIGNED_FORMAT, NetworkitBinaryWeights::SIGNED_FORMAT,
          NetworkitBinaryWeights::DOUBLE_FORMAT, NetworkitBinaryWeights::FLOAT_FORMAT}) {
        const GraphT G = TestFixture::smallGraph(true, false);
        const std::vector<uint8_t> buffer =
            GenericNetworkitBinaryWriter<GraphT>(32, format).writeToBuffer(G);
        TestFixture::expectSameGraph(G,
                                     GenericNetworkitBinaryReader<GraphT>{}.readFromBuffer(buffer));
    }
}

// The upper edge id bound is one more than the largest edge id, so that a new edge gets an unused
// id.
TYPED_TEST(NetworkitBinaryIOGTest, testReadEdgeIdsAllowNewEdges) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    GraphT G = TestFixture::smallGraph(false, false);
    G.indexEdges();
    const std::vector<uint8_t> buffer = GenericNetworkitBinaryWriter<GraphT>{}.writeToBuffer(G);
    GraphT read = GenericNetworkitBinaryReader<GraphT>{}.readFromBuffer(buffer);

    read.addEdge(NodeT{2}, NodeT{3});
    const edgeid newId = read.edgeId(NodeT{2}, NodeT{3});
    count edgesWithNewId = 0;
    read.forEdges([&](NodeT, NodeT, EdgeWeightT, edgeid id) -> void {
        if (id == newId) {
            ++edgesWithNewId;
        }
    });
    EXPECT_EQ(edgesWithNewId, count{1});
}

TYPED_TEST(NetworkitBinaryIOGTest, testReadWeightsThatDoNotFitIntegralTypes) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const auto fileWithWeight = [](double weight) -> std::vector<uint8_t> {
        Graph G(2, true);
        G.addEdge(0, 1, weight);
        return NetworkitBinaryWriter{}.writeToBuffer(G);
    };
    GenericNetworkitBinaryReader<GraphT> reader;

    // 0.5 is written as float, -4 as signed varint, and 256 as unsigned varint.
    if constexpr (std::is_integral_v<EdgeWeightT>) {
        EXPECT_THAT([&]() -> void { reader.readFromBuffer(fileWithWeight(0.5)); },
                    ThrowsMessage<std::runtime_error>(HasSubstr("cannot represent")));
    } else {
        EXPECT_EQ(reader.readFromBuffer(fileWithWeight(0.5)).weight(NodeT{0}, NodeT{1}),
                  EdgeWeightT{0.5});
    }
    if constexpr (std::is_unsigned_v<EdgeWeightT>) {
        EXPECT_THROW(reader.readFromBuffer(fileWithWeight(-4.0)), std::runtime_error);
    } else {
        EXPECT_EQ(reader.readFromBuffer(fileWithWeight(-4.0)).weight(NodeT{0}, NodeT{1}),
                  EdgeWeightT{-4});
    }
    if constexpr (std::is_same_v<EdgeWeightT, uint8_t>) {
        EXPECT_THROW(reader.readFromBuffer(fileWithWeight(256.0)), std::runtime_error);
    } else {
        EXPECT_EQ(reader.readFromBuffer(fileWithWeight(256.0)).weight(NodeT{0}, NodeT{1}),
                  EdgeWeightT{256});
    }
}

TEST(NetworkitBinaryIOGTest, testReadLargeIntegerWeightsExactly) {
    using GraphT = AdjListGraph<uint32_t, uint64_t>;
    // Above 2^53, so that a detour through double would round it.
    constexpr uint64_t weight = (uint64_t{1} << 60) + 1;
    GraphT G(2, true);
    G.addEdge(0u, 1u, weight);
    const std::vector<uint8_t> buffer = GenericNetworkitBinaryWriter<GraphT>{}.writeToBuffer(G);
    EXPECT_EQ(GenericNetworkitBinaryReader<GraphT>{}.readFromBuffer(buffer).weight(0u, 1u), weight);
}

TEST(NetworkitBinaryIOGTest, testReadThrowsIfTooManyNodes) {
    using Reader = GenericNetworkitBinaryReader<AdjListGraph<uint8_t, uint8_t>>;
    using Writer = GenericNetworkitBinaryWriter<AdjListGraph<uint16_t, uint8_t>>;
    constexpr count maxNodes = std::numeric_limits<uint8_t>::max();

    const std::vector<uint8_t> fits =
        Writer{}.writeToBuffer(AdjListGraph<uint16_t, uint8_t>(maxNodes));
    EXPECT_EQ(Reader{}.readFromBuffer(fits).numberOfNodes(), maxNodes);

    const std::vector<uint8_t> tooMany =
        Writer{}.writeToBuffer(AdjListGraph<uint16_t, uint8_t>(maxNodes + 1));
    EXPECT_THAT([&]() -> void { Reader{}.readFromBuffer(tooMany); },
                ThrowsMessage<std::runtime_error>(HasSubstr("can hold at most")));
}

} // namespace
} // namespace NetworKit

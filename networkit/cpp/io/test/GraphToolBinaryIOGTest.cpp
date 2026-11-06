#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/GraphToolBinaryReader.hpp>
#include <networkit/io/GraphToolBinaryWriter.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

template <class GraphT_>
class GraphToolBinaryIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
};

using GraphToolBinaryIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(GraphToolBinaryIOGTest, GraphToolBinaryIOTestTypes,
                 /*Comma needed for variadic macro.*/);

TYPED_TEST(GraphToolBinaryIOGTest, testWriteAndReadBack) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;

    for (const bool directed : {false, true}) {
        for (const bool littleEndian : {false, true}) {
            GraphT G(13, false, directed);
            G.addEdge(NodeT{0}, NodeT{1});
            G.addEdge(NodeT{1}, NodeT{12});
            G.addEdge(NodeT{12}, NodeT{11});
            G.addEdge(NodeT{11}, NodeT{0});
            G.addEdge(NodeT{11}, NodeT{11});
            // The file has continuous node ids, so the ids after the deleted node shift by one.
            G.removeNode(NodeT{5});
            const auto readId = [](NodeT u) -> NodeT {
                return u < NodeT{5} ? u : static_cast<NodeT>(u - 1);
            };

            const std::string path = directory.path("graph.gt");
            GenericGraphToolBinaryWriter<GraphT>(littleEndian).write(G, path);
            const GraphT read = GenericGraphToolBinaryReader<GraphT>{}.read(path);
            EXPECT_EQ(read.isDirected(), directed);
            EXPECT_FALSE(read.isWeighted());
            EXPECT_TRUE(read.checkConsistency());
            EXPECT_EQ(read.numberOfNodes(), count{12});
            EXPECT_EQ(read.numberOfEdges(), G.numberOfEdges());
            G.forEdges(
                [&](NodeT u, NodeT v) -> void { EXPECT_TRUE(read.hasEdge(readId(u), readId(v))); });
        }
    }
}

TEST(GraphToolBinaryIOGTest, testReadThrowsIfTooManyNodes) {
    constexpr count maxNodes = std::numeric_limits<uint8_t>::max();
    using Reader = GenericGraphToolBinaryReader<AdjListGraph<uint8_t, uint8_t>>;
    using Writer = GenericGraphToolBinaryWriter<AdjListGraph<uint16_t, uint8_t>>;
    const TemporaryDirectory directory;

    const std::string fits = directory.path("fits.gt");
    Writer{}.write(AdjListGraph<uint16_t, uint8_t>(maxNodes), fits);
    EXPECT_EQ(Reader{}.read(fits).numberOfNodes(), maxNodes);

    const std::string tooMany = directory.path("too_many.gt");
    Writer{}.write(AdjListGraph<uint16_t, uint8_t>(maxNodes + 1), tooMany);
    EXPECT_THAT([&]() -> void { Reader{}.read(tooMany); },
                ThrowsMessage<std::runtime_error>(HasSubstr("can hold at most")));
}

TEST(GraphToolBinaryIOGTest, testReadThrowsIfNodeIdIsTooLarge) {
    // Header, empty comment, undirected, two nodes; node 0 has the neighbor 5.
    std::string file("\xe2\x9b\xbe\x20\x67\x74\x01\x00", 8);
    file += std::string(8, '\0') + '\0';
    file += std::string("\x02", 1) + std::string(7, '\0');
    file += std::string("\x01", 1) + std::string(7, '\0') + '\x05';
    file += std::string(8, '\0');

    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("large_id.gt", file);
    EXPECT_THAT([&]() -> void { GraphToolBinaryReader{}.read(path); },
                ThrowsMessage<std::runtime_error>(HasSubstr("valid ids are 0 to 1")));
}

} // namespace
} // namespace NetworKit

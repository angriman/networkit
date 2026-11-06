#ifndef NETWORKIT_WINDOWS

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <netinet/in.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/DibapGraphReader.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

// Returns a DibaP GI file of an unweighted graph with V nodes, adjacency offsets of and targets to,
// dvw vertex weights per node, and dxy coordinates per node.
std::string dibapFile(int V, int dvw, const std::vector<int> &of, const std::vector<int> &to,
                      int dxy, const std::vector<float> &xy) {
    std::string file;
    const auto putShort = [&](uint16_t value) -> void {
        const uint16_t bigEndian = htons(value);
        file.append(reinterpret_cast<const char *>(&bigEndian), sizeof(bigEndian));
    };
    const auto putInt = [&](int value) -> void {
        const uint32_t bigEndian = htonl(static_cast<uint32_t>(value));
        file.append(reinterpret_cast<const char *>(&bigEndian), sizeof(bigEndian));
    };
    putShort(('G' << 8) | 'I');
    putInt(V);
    putInt(dvw);
    for (int i = 0; i < V * dvw; ++i) {
        putInt(1);
    }
    for (const int offset : of) {
        putInt(offset);
    }
    for (const int target : to) {
        putInt(target);
    }
    putInt(0); // no edge weights
    putInt(dxy);
    file.append(reinterpret_cast<const char *>(xy.data()), xy.size() * sizeof(float));
    return file;
}

template <class GraphT_>
class DibapGraphIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
};

using DibapGraphIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(DibapGraphIOGTest, DibapGraphIOTestTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(DibapGraphIOGTest, testRead) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    // The path 0 - 1 - 2, with two vertex weights and two coordinates per node.
    const std::string path = directory.writeFile(
        "path.gi", dibapFile(3, 2, {0, 1, 3, 4}, {1, 0, 2, 1}, 2, {0.f, 1.f, 2.f, 3.f, 4.f, 5.f}));
    GenericDibapGraphReader<GraphT> reader;
    const GraphT G = reader.read(path);
    EXPECT_FALSE(G.isDirected());
    EXPECT_FALSE(G.isWeighted());
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{3});
    EXPECT_EQ(G.numberOfEdges(), count{2});
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{1}));
    EXPECT_TRUE(G.hasEdge(NodeT{1}, NodeT{2}));

    const auto &coordinates = reader.getCoordinates();
    ASSERT_EQ(coordinates.size(), size_t{3});
    EXPECT_EQ(coordinates[2][0], 4.f);
    EXPECT_EQ(coordinates[2][1], 5.f);
}

TYPED_TEST(DibapGraphIOGTest, testReadAirfoil) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const Graph expected = DibapGraphReader{}.read("input/airfoil1.gi");
    GenericDibapGraphReader<GraphT> reader;
    if (expected.numberOfNodes() > static_cast<count>(std::numeric_limits<NodeT>::max())) {
        EXPECT_THAT([&]() -> void { reader.read("input/airfoil1.gi"); },
                    ThrowsMessage<std::runtime_error>(HasSubstr("can hold at most")));
        return;
    }
    const GraphT G = reader.read("input/airfoil1.gi");
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), expected.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), expected.numberOfEdges());
    expected.forEdges([&](node u, node v) -> void {
        EXPECT_TRUE(G.hasEdge(static_cast<NodeT>(u), static_cast<NodeT>(v)));
    });
}

TEST(DibapGraphIOGTest, testReadThrowsIfNodeIdIsOutOfRange) {
    const TemporaryDirectory directory;
    const std::string path =
        directory.writeFile("large_id.gi", dibapFile(2, 0, {0, 1, 2}, {1, 2}, 0, {}));
    EXPECT_THAT([&]() -> void { DibapGraphReader{}.read(path); },
                ThrowsMessage<std::runtime_error>(HasSubstr("valid ids are 0 to 1")));
}

} // namespace
} // namespace NetworKit

#endif // NETWORKIT_WINDOWS

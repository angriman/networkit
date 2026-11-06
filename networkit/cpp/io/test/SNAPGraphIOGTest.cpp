#include <cstdint>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/io/SNAPGraphWriter.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

template <class GraphT_>
class SNAPGraphIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

protected:
    static std::string readFile(const std::string &path) {
        std::ifstream file(path);
        return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
    }
};

using SNAPGraphIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(SNAPGraphIOGTest, SNAPGraphIOTestTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(SNAPGraphIOGTest, testReadRemapped) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("remapped.txt", "# A comment\n"
                                                                 "20 10\n"
                                                                 "\n"
                                                                 "10\t30\r\n"
                                                                 "30 20\n");
    const GraphT G = GenericSNAPGraphReader<GraphT>{}.read(path);
    EXPECT_FALSE(G.isDirected());
    EXPECT_FALSE(G.isWeighted());
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{3});
    EXPECT_EQ(G.numberOfEdges(), count{3});
    // The ids are assigned in order of appearance: 20, 10, 30.
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{1}));
    EXPECT_TRUE(G.hasEdge(NodeT{1}, NodeT{2}));
    EXPECT_TRUE(G.hasEdge(NodeT{2}, NodeT{0}));
}

TYPED_TEST(SNAPGraphIOGTest, testReadDirectedWithoutRemapping) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("directed.txt", "0 1\n"
                                                                 "1 0\n"
                                                                 "3 1\n");
    const GraphT G = GenericSNAPGraphReader<GraphT>(true, false).read(path);
    EXPECT_TRUE(G.isDirected());
    EXPECT_TRUE(G.checkConsistency());
    // Node 2 is isolated.
    EXPECT_EQ(G.numberOfNodes(), count{4});
    EXPECT_EQ(G.numberOfEdges(), count{3});
    EXPECT_TRUE(G.hasEdge(NodeT{3}, NodeT{1}));
    EXPECT_FALSE(G.hasEdge(NodeT{1}, NodeT{3}));
}

TYPED_TEST(SNAPGraphIOGTest, testReadTwiceWithSameReader) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    GenericSNAPGraphReader<GraphT> reader;
    EXPECT_EQ(reader.read(directory.writeFile("first.txt", "5 6\n6 7\n")).numberOfNodes(),
              count{3});

    // Node 7 is new in this file, so it must get id 1 and not the id 2 of the first file.
    const GraphT G = reader.read(directory.writeFile("second.txt", "6 7\n"));
    EXPECT_EQ(G.numberOfNodes(), count{2});
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{1}));
}

TYPED_TEST(SNAPGraphIOGTest, testReadThrowsIfLineIsMalformed) {
    using GraphT = typename TestFixture::GraphT;
    const TemporaryDirectory directory;
    GenericSNAPGraphReader<GraphT> reader;
    EXPECT_THROW(reader.read(directory.writeFile("letters.txt", "0 x\n")), std::runtime_error);
    EXPECT_THROW(reader.read(directory.writeFile("negative.txt", "0 -1\n")), std::runtime_error);
    EXPECT_THROW(reader.read(directory.writeFile("fraction.txt", "0 1.5\n")), std::runtime_error);
    EXPECT_THROW(reader.read(directory.writeFile("single.txt", "0\n")), std::runtime_error);
}

TEST(SNAPGraphIOGTest, testReadThrowsIfNodeIdDoesNotFitNodeType) {
    using Reader = GenericSNAPGraphReader<AdjListGraph<uint16_t, double>>;
    constexpr uint32_t maxNode = std::numeric_limits<uint16_t>::max();
    const TemporaryDirectory directory;

    // The largest node id maxNode - 1 is still smaller than NullNodeId<uint16_t>.
    const std::string fits =
        directory.writeFile("max_id.txt", "0 " + std::to_string(maxNode - 1) + "\n");
    EXPECT_EQ(Reader(false, false).read(fits).numberOfNodes(), count{maxNode});

    const std::string nullId =
        directory.writeFile("null_id.txt", "0 " + std::to_string(maxNode) + "\n");
    EXPECT_THAT([&]() -> void { Reader(false, false).read(nullId); },
                ThrowsMessage<std::runtime_error>(HasSubstr("does not fit")));

    // With remapping, only the number of nodes matters.
    EXPECT_EQ(Reader{}.read(nullId).numberOfNodes(), count{2});
}

TEST(SNAPGraphIOGTest, testReadThrowsIfTooManyNodes) {
    using Reader = GenericSNAPGraphReader<AdjListGraph<uint8_t, uint8_t>>;
    // A path through the largest number of nodes that fits, 255, and then one more.
    constexpr uint32_t maxNodes = std::numeric_limits<uint8_t>::max();
    std::string path;
    for (uint32_t u = 0; u + 1 < maxNodes; ++u) {
        path += std::to_string(1000 + u) + " " + std::to_string(1000 + u + 1) + "\n";
    }
    const TemporaryDirectory directory;
    EXPECT_EQ(Reader{}.read(directory.writeFile("fits.txt", path)).numberOfNodes(),
              count{maxNodes});

    path += std::to_string(1000 + maxNodes - 1) + " " + std::to_string(1000 + maxNodes) + "\n";
    EXPECT_THAT([&]() -> void { Reader{}.read(directory.writeFile("too_many.txt", path)); },
                ThrowsMessage<std::runtime_error>(HasSubstr("more nodes")));
}

TYPED_TEST(SNAPGraphIOGTest, testWrite) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;

    // Node ids from 10 on are printed as text also for 8-bit node types.
    GraphT G(12);
    G.addEdge(NodeT{0}, NodeT{1});
    G.addEdge(NodeT{1}, NodeT{11});
    G.addEdge(NodeT{10}, NodeT{11});

    const std::string path = directory.path("graph.gr");
    GenericSNAPGraphWriter<GraphT>{}.write(G, path);
    EXPECT_EQ(TestFixture::readFile(path), "p 12 3 u u 0\n"
                                           "1 0\n"
                                           "11 1\n"
                                           "11 10\n");
}

} // namespace
} // namespace NetworKit

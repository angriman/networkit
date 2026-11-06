#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/MTXGraphReader.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

template <class GraphT_>
class MTXGraphIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
};

using MTXGraphIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(MTXGraphIOGTest, MTXGraphIOTestTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(MTXGraphIOGTest, testReadPatternSymmetric) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("pattern.mtx", "%%MatrixMarket matrix coordinate "
                                                                "pattern symmetric\n"
                                                                "% A comment\n"
                                                                "3 3 2\n"
                                                                "2 1\n"
                                                                "\n"
                                                                "3 1\r\n");
    const GraphT G = GenericMTXGraphReader<GraphT>{}.read(path);
    EXPECT_FALSE(G.isDirected());
    EXPECT_FALSE(G.isWeighted());
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{3});
    EXPECT_EQ(G.numberOfEdges(), count{2});
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{1}));
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{2}));
}

TYPED_TEST(MTXGraphIOGTest, testReadIntegerGeneral) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    // A rectangular matrix: the graph has a node per row or column, whichever are more.
    const std::string path = directory.writeFile("integer.mtx", "%%MatrixMarket matrix coordinate "
                                                                "integer general\n"
                                                                "2 4 3\n"
                                                                "1 2 7\n"
                                                                "2 1 100\n"
                                                                "2 4 3\n");
    const GraphT G = GenericMTXGraphReader<GraphT>{}.read(path);
    EXPECT_TRUE(G.isDirected());
    EXPECT_TRUE(G.isWeighted());
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{4});
    EXPECT_EQ(G.numberOfEdges(), count{3});
    EXPECT_EQ(G.weight(NodeT{0}, NodeT{1}), EdgeWeightT{7});
    EXPECT_EQ(G.weight(NodeT{1}, NodeT{0}), EdgeWeightT{100});
    EXPECT_EQ(G.weight(NodeT{1}, NodeT{3}), EdgeWeightT{3});
    EXPECT_FALSE(G.hasEdge(NodeT{3}, NodeT{1}));
}

TYPED_TEST(MTXGraphIOGTest, testReadFractionalWeights) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    // Weights exactly representable by float and double.
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("real.mtx", "%%MatrixMarket matrix coordinate "
                                                             "real symmetric\n"
                                                             "2 2 1\n"
                                                             "2 1 0.5\n");
    GenericMTXGraphReader<GraphT> reader;
    if constexpr (std::is_integral_v<EdgeWeightT>) {
        EXPECT_THAT([&]() -> void { reader.read(path); },
                    ThrowsMessage<std::runtime_error>(HasSubstr("not valid")));
    } else {
        EXPECT_EQ(reader.read(path).weight(NodeT{0}, NodeT{1}), EdgeWeightT{0.5});
    }
}

TYPED_TEST(MTXGraphIOGTest, testReadThrowsIfNodeIdIsOutOfRange) {
    using GraphT = typename TestFixture::GraphT;
    const TemporaryDirectory directory;
    const std::string header = "%%MatrixMarket matrix coordinate pattern general\n"
                               "2 2 1\n";
    GenericMTXGraphReader<GraphT> reader;
    EXPECT_THAT([&]() -> void { reader.read(directory.writeFile("zero.mtx", header + "0 1\n")); },
                ThrowsMessage<std::runtime_error>(HasSubstr("start at 1")));
    EXPECT_THAT([&]() -> void { reader.read(directory.writeFile("large.mtx", header + "1 3\n")); },
                ThrowsMessage<std::runtime_error>(HasSubstr("only 2 rows")));
    EXPECT_THROW(reader.read(directory.writeFile("negative.mtx", header + "1 -1\n")),
                 std::runtime_error);
    EXPECT_THROW(reader.read(directory.writeFile("single.mtx", header + "1\n")),
                 std::runtime_error);
}

TEST(MTXGraphIOGTest, testReadThrowsIfMissingWeight) {
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("missing.mtx", "%%MatrixMarket matrix coordinate "
                                                                "real general\n"
                                                                "2 2 1\n"
                                                                "1 2\n");
    EXPECT_THAT([&]() -> void { MTXGraphReader{}.read(path); },
                ThrowsMessage<std::runtime_error>(HasSubstr("weight was expected")));
}

TEST(MTXGraphIOGTest, testReadLargeIntegerWeightsExactly) {
    using GraphT = AdjListGraph<uint32_t, uint64_t>;
    // 2^53 + 1 is not representable by double.
    constexpr uint64_t weight = (uint64_t{1} << 53) + 1;
    const TemporaryDirectory directory;
    const std::string path =
        directory.writeFile("large.mtx", "%%MatrixMarket matrix coordinate integer general\n"
                                         "2 2 1\n"
                                         "1 2 "
                                             + std::to_string(weight) + "\n");
    EXPECT_EQ(GenericMTXGraphReader<GraphT>{}.read(path).weight(0u, 1u), weight);
}

TEST(MTXGraphIOGTest, testReadThrowsIfTooManyNodes) {
    using Reader = GenericMTXGraphReader<AdjListGraph<uint8_t, uint8_t>>;
    constexpr count maxNodes = std::numeric_limits<uint8_t>::max();
    const auto mtx = [](count n) -> std::string {
        return "%%MatrixMarket matrix coordinate pattern general\n" + std::to_string(n) + " "
               + std::to_string(n) + " 0\n";
    };
    const TemporaryDirectory directory;
    EXPECT_EQ(Reader{}.read(directory.writeFile("fits.mtx", mtx(maxNodes))).numberOfNodes(),
              maxNodes);
    EXPECT_THAT(
        [&]() -> void { Reader{}.read(directory.writeFile("many.mtx", mtx(maxNodes + 1))); },
        ThrowsMessage<std::runtime_error>(HasSubstr("can hold at most")));
}

} // namespace
} // namespace NetworKit

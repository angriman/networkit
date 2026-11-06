#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/KONECTGraphReader.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

template <class GraphT_>
class KONECTGraphIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
};

using KONECTGraphIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(KONECTGraphIOGTest, KONECTGraphIOTestTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(KONECTGraphIOGTest, testReadUnweighted) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("unweighted.konect", "% sym unweighted\n"
                                                                      "1 2\n"
                                                                      "% A comment\n"
                                                                      "\n"
                                                                      "3\t2\r\n");
    const GraphT G = GenericKONECTGraphReader<GraphT>{}.read(path);
    EXPECT_FALSE(G.isDirected());
    EXPECT_FALSE(G.isWeighted());
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{3});
    EXPECT_EQ(G.numberOfEdges(), count{2});
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{1}));
    EXPECT_TRUE(G.hasEdge(NodeT{1}, NodeT{2}));
}

TYPED_TEST(KONECTGraphIOGTest, testReadWeightedWithNumberOfNodes) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    // Node 4 is isolated.
    const std::string path = directory.writeFile("weighted.konect", "% asym weighted\n"
                                                                    "% 3 4\n"
                                                                    "1 2 7\n"
                                                                    "2 1 100\n"
                                                                    "3 2 3\n");
    const GraphT G = GenericKONECTGraphReader<GraphT>{}.read(path);
    EXPECT_TRUE(G.isDirected());
    EXPECT_TRUE(G.isWeighted());
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{4});
    EXPECT_EQ(G.numberOfEdges(), count{3});
    EXPECT_EQ(G.weight(NodeT{0}, NodeT{1}), EdgeWeightT{7});
    EXPECT_EQ(G.weight(NodeT{1}, NodeT{0}), EdgeWeightT{100});
    EXPECT_EQ(G.weight(NodeT{2}, NodeT{1}), EdgeWeightT{3});
}

TYPED_TEST(KONECTGraphIOGTest, testReadMultipleEdges) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using Reader = GenericKONECTGraphReader<GraphT>;
    const TemporaryDirectory directory;

    const std::string positive = directory.writeFile("positive.konect", "% sym positive\n"
                                                                        "1 2\n"
                                                                        "1 2\n"
                                                                        "2 3\n");
    const GraphT discarded = Reader(false, GraphReaderOptions::DISCARD_EDGES).read(positive);
    EXPECT_FALSE(discarded.isWeighted());
    EXPECT_EQ(discarded.numberOfEdges(), count{2});
    const GraphT summed = Reader(false, GraphReaderOptions::SUM_WEIGHTS_UP).read(positive);
    EXPECT_TRUE(summed.isWeighted());
    EXPECT_EQ(summed.weight(NodeT{0}, NodeT{1}), EdgeWeightT{2});
    EXPECT_EQ(summed.weight(NodeT{1}, NodeT{2}), EdgeWeightT{1});

    const std::string weighted = directory.writeFile("multiweighted.konect", "% sym multiweighted\n"
                                                                             "1 2 5\n"
                                                                             "1 2 3\n"
                                                                             "1 2 4\n");
    const GraphT minimum = Reader(false, GraphReaderOptions::KEEP_MINIMUM_WEIGHT).read(weighted);
    EXPECT_EQ(minimum.numberOfEdges(), count{1});
    EXPECT_EQ(minimum.weight(NodeT{0}, NodeT{1}), EdgeWeightT{3});
}

TYPED_TEST(KONECTGraphIOGTest, testReadFractionalWeights) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    // Weights exactly representable by float and double.
    const std::string path = directory.writeFile("fractional.konect", "% sym weighted\n"
                                                                      "1 2 0.5\n");
    GenericKONECTGraphReader<GraphT> reader;
    if constexpr (std::is_integral_v<EdgeWeightT>) {
        EXPECT_THAT([&]() -> void { reader.read(path); },
                    ThrowsMessage<std::runtime_error>(HasSubstr("Scanning weight 0.5 failed")));
    } else {
        EXPECT_EQ(reader.read(path).weight(NodeT{0}, NodeT{1}), EdgeWeightT{0.5});
    }
}

TYPED_TEST(KONECTGraphIOGTest, testReadRemapped) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("remapped.konect", "% asym unweighted\n"
                                                                    "20 10\n"
                                                                    "10 30\n");
    const GraphT G = GenericKONECTGraphReader<GraphT>(true).read(path);
    EXPECT_EQ(G.numberOfNodes(), count{3});
    // The ids are assigned in order of appearance: 20, 10, 30.
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{1}));
    EXPECT_TRUE(G.hasEdge(NodeT{1}, NodeT{2}));
}

// The file claims 2 nodes but has 4: the reader falls back to remapping, and must keep the ids of
// the nodes it has read so far.
TYPED_TEST(KONECTGraphIOGTest, testReadWrongNumberOfNodesWithRemapping) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("wrong_count.konect", "% asym unweighted\n"
                                                                       "% 3 2\n"
                                                                       "1 2\n"
                                                                       "2 5\n"
                                                                       "2 3\n");
    const GraphT G = GenericKONECTGraphReader<GraphT>(true).read(path);
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{4});
    EXPECT_EQ(G.numberOfEdges(), count{3});
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{1}));
    EXPECT_TRUE(G.hasEdge(NodeT{1}, NodeT{2}));
    EXPECT_TRUE(G.hasEdge(NodeT{1}, NodeT{3}));
}

TYPED_TEST(KONECTGraphIOGTest, testReadThrowsIfNodeIdIsInvalid) {
    using GraphT = typename TestFixture::GraphT;
    const TemporaryDirectory directory;
    GenericKONECTGraphReader<GraphT> reader;
    EXPECT_THAT(
        [&]() -> void {
            reader.read(directory.writeFile("zero.konect", "% sym unweighted\n0 1\n"));
        },
        ThrowsMessage<std::runtime_error>(HasSubstr("start at 1")));
    EXPECT_THROW(reader.read(directory.writeFile("negative.konect", "% sym unweighted\n1 -2\n")),
                 std::runtime_error);
}

TEST(KONECTGraphIOGTest, testReadThrowsIfNodeIdDoesNotFitNodeType) {
    using Reader = GenericKONECTGraphReader<AdjListGraph<uint16_t, double>>;
    constexpr uint32_t maxNode = std::numeric_limits<uint16_t>::max();
    const TemporaryDirectory directory;

    // The largest 1-based node id maxNode becomes maxNode - 1, which is smaller than
    // NullNodeId<uint16_t>.
    const std::string fits = directory.writeFile(
        "max_id.konect", "% sym unweighted\n1 " + std::to_string(maxNode) + "\n");
    EXPECT_EQ(Reader{}.read(fits).numberOfNodes(), count{maxNode});

    const std::string tooLarge = directory.writeFile(
        "large_id.konect", "% sym unweighted\n1 " + std::to_string(maxNode + 1) + "\n");
    EXPECT_THAT([&]() -> void { Reader{}.read(tooLarge); },
                ThrowsMessage<std::runtime_error>(HasSubstr("does not fit")));

    // With remapping, only the number of nodes matters.
    EXPECT_EQ(Reader(true).read(tooLarge).numberOfNodes(), count{2});
}

TEST(KONECTGraphIOGTest, testReadThrowsIfTooManyNodes) {
    using Reader = GenericKONECTGraphReader<AdjListGraph<uint8_t, uint8_t>>;
    constexpr uint32_t maxNodes = std::numeric_limits<uint8_t>::max();
    const TemporaryDirectory directory;

    const auto withNumberOfNodes = [](uint32_t n) -> std::string {
        return "% sym unweighted\n% 1 " + std::to_string(n) + "\n1 2\n";
    };
    EXPECT_EQ(Reader{}
                  .read(directory.writeFile("fits.konect", withNumberOfNodes(maxNodes)))
                  .numberOfNodes(),
              count{maxNodes});
    EXPECT_THAT(
        [&]() -> void {
            Reader{}.read(directory.writeFile("many.konect", withNumberOfNodes(maxNodes + 1)));
        },
        ThrowsMessage<std::runtime_error>(HasSubstr("can hold at most")));

    // With remapping, a path through 256 distinct ids.
    std::string path = "% sym unweighted\n";
    for (uint32_t u = 0; u < maxNodes; ++u) {
        path += std::to_string(1000 + u) + " " + std::to_string(1000 + u + 1) + "\n";
    }
    EXPECT_THAT([&]() -> void { Reader(true).read(directory.writeFile("remap.konect", path)); },
                ThrowsMessage<std::runtime_error>(HasSubstr("more nodes")));
}

} // namespace
} // namespace NetworKit

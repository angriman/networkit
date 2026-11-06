#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/RBGraphReader.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

// Returns an RB file of an n x n matrix in compressed column format, with one line each for the
// pointers, the row indices, and (unless empty) the values.
std::string rbFile(const std::string &type, count n, const std::string &pointers,
                   const std::string &rows, count nonzeros, const std::string &values) {
    const count valueLines = values.empty() ? 0 : 1;
    std::string text =
        "A test matrix                                                           KEY\n";
    text += std::to_string(2 + valueLines) + " 1 1 " + std::to_string(valueLines) + " 0\n";
    text += type + " " + std::to_string(n) + " " + std::to_string(n) + " "
            + std::to_string(nonzeros) + "\n";
    text += "(10I8) (10I8) (5E15.8)\n";
    text += pointers + "\n" + rows + "\n";
    if (!values.empty()) {
        text += values + "\n";
    }
    return text;
}

template <class GraphT_>
class RBGraphIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
};

using RBGraphIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint8_t, uint8_t>>;

TYPED_TEST_SUITE(RBGraphIOGTest, RBGraphIOTestTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(RBGraphIOGTest, testReadPatternSymmetric) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    // The entries (2, 1) and (3, 1) of the lower triangle.
    const std::string path =
        directory.writeFile("pattern.rb", rbFile("PSA", 3, "1 3 3 3", "2 3", 2, ""));
    const GraphT G = GenericRBGraphReader<GraphT>{}.read(path);
    EXPECT_FALSE(G.isDirected());
    EXPECT_FALSE(G.isWeighted());
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{3});
    EXPECT_EQ(G.numberOfEdges(), count{2});
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{1}));
    EXPECT_TRUE(G.hasEdge(NodeT{0}, NodeT{2}));
}

TYPED_TEST(RBGraphIOGTest, testReadIntegerUnsymmetric) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    // The entries (2, 1) = 7, (1, 2) = 100, and (2, 2) = 3.
    const std::string path =
        directory.writeFile("integer.rb", rbFile("IUA", 2, "1 2 4", "2 1 2", 3, "7 100 3"));
    const GraphT G = GenericRBGraphReader<GraphT>{}.read(path);
    EXPECT_TRUE(G.isDirected());
    EXPECT_TRUE(G.isWeighted());
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(G.numberOfNodes(), count{2});
    EXPECT_EQ(G.numberOfEdges(), count{3});
    EXPECT_EQ(G.weight(NodeT{1}, NodeT{0}), EdgeWeightT{7});
    EXPECT_EQ(G.weight(NodeT{0}, NodeT{1}), EdgeWeightT{100});
    EXPECT_EQ(G.weight(NodeT{1}, NodeT{1}), EdgeWeightT{3});
}

TYPED_TEST(RBGraphIOGTest, testReadWeightsThatDoNotFitIntegralTypes) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    GenericRBGraphReader<GraphT> reader;
    // Weights exactly representable by float and double.
    const std::string fraction =
        directory.writeFile("fraction.rb", rbFile("RSA", 2, "1 2 2", "2", 1, "0.5"));
    const std::string negative =
        directory.writeFile("negative.rb", rbFile("RSA", 2, "1 2 2", "2", 1, "-4.0"));
    const std::string large =
        directory.writeFile("large.rb", rbFile("RSA", 2, "1 2 2", "2", 1, "256.0"));

    if constexpr (std::is_integral_v<EdgeWeightT>) {
        EXPECT_THAT([&]() -> void { reader.read(fraction); },
                    ThrowsMessage<std::runtime_error>(HasSubstr("not valid")));
    } else {
        EXPECT_EQ(reader.read(fraction).weight(NodeT{0}, NodeT{1}), EdgeWeightT{0.5});
    }
    if constexpr (std::is_unsigned_v<EdgeWeightT>) {
        EXPECT_THROW(reader.read(negative), std::runtime_error);
    } else {
        EXPECT_EQ(reader.read(negative).weight(NodeT{0}, NodeT{1}), EdgeWeightT{-4});
    }
    if constexpr (std::is_same_v<EdgeWeightT, uint8_t>) {
        EXPECT_THROW(reader.read(large), std::runtime_error);
    } else {
        EXPECT_EQ(reader.read(large).weight(NodeT{0}, NodeT{1}), EdgeWeightT{256});
    }
}

TYPED_TEST(RBGraphIOGTest, testReadThrowsIfRowIndexIsOutOfRange) {
    using GraphT = typename TestFixture::GraphT;
    const TemporaryDirectory directory;
    const std::string path =
        directory.writeFile("large_row.rb", rbFile("PUA", 2, "1 2 2", "3", 1, ""));
    EXPECT_THAT([&]() -> void { GenericRBGraphReader<GraphT>{}.read(path); },
                ThrowsMessage<std::runtime_error>(HasSubstr("valid ones are 1 to 2")));
}

TEST(RBGraphIOGTest, testReadThrowsIfTooManyNodes) {
    using Reader = GenericRBGraphReader<AdjListGraph<uint8_t, uint8_t>>;
    constexpr count maxNodes = std::numeric_limits<uint8_t>::max();
    // A matrix without entries: all pointers are 1.
    const auto emptyMatrix = [](count n) -> std::string {
        std::string pointers = "1";
        for (count i = 0; i < n; ++i) {
            pointers += " 1";
        }
        return rbFile("PUA", n, pointers, "", 0, "");
    };
    const TemporaryDirectory directory;
    EXPECT_EQ(Reader{}.read(directory.writeFile("fits.rb", emptyMatrix(maxNodes))).numberOfNodes(),
              maxNodes);
    EXPECT_THAT(
        [&]() -> void { Reader{}.read(directory.writeFile("many.rb", emptyMatrix(maxNodes + 1))); },
        ThrowsMessage<std::runtime_error>(HasSubstr("can hold at most")));
}

} // namespace
} // namespace NetworKit

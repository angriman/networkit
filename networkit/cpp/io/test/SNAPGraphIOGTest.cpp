#include <cstdint>
#include <fstream>
#include <iterator>
#include <string>

#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/SNAPGraphWriter.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

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

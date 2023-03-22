#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <networkit/graph/CSRGraph.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

struct GraphProperties {
    bool isDirected;
};

class CSRWeightedGraphGTest : public ::testing::TestWithParam<GraphProperties> {
protected:
    bool isDirected() const noexcept { return GetParam().isDirected; }
};

TEST_P(CSRWeightedGraphGTest, DefaultCtorCreatesEmptyGraph) {
    CSRGraph<int, float> G(0, isDirected());

    EXPECT_EQ(G.numberOfNodes(), 0);
    EXPECT_EQ(G.numberOfEdges(), 0);
}

INSTANTIATE_TEST_SUITE_P(TestCSRWeightedGraph, CSRWeightedGraphGTest,
                         testing::Values(GraphProperties{false}, GraphProperties{true}));

} // namespace NetworKit

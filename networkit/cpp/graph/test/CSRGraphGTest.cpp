#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <networkit/graph/CSRGraph.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace {

using ::testing::Each;
using ::testing::TestWithParam;

struct GraphProperties {
    bool isDirected;
};

class CSRGraphGTest : public TestWithParam<GraphProperties> {
protected:
    bool isDirected() const noexcept { return GetParam().isDirected; }
};

TEST_P(CSRGraphGTest, DefaultCtorCreatesEmptyGraph) {
    CSRGraph<int> G(0, isDirected());

    EXPECT_EQ(G.numberOfNodes(), 0);
    EXPECT_EQ(G.numberOfEdges(), 0);
}

TEST_P(CSRGraphGTest, DefaultCtorCreatesGraphWithCorrectSize) {
    CSRGraph<int> G(5, isDirected());

    EXPECT_EQ(G.numberOfNodes(), 5);
    EXPECT_EQ(G.numberOfEdges(), 0);
}

TEST_P(CSRGraphGTest, CtorFromDynGraphCreatesGraphWithCorrectSize) {
    Graph dynG(10, false, isDirected());
    dynG.addEdge(0, 1);
    dynG.addEdge(1, 2);
    dynG.addEdge(2, 3);

    CSRGraph<int> G(dynG);
    EXPECT_EQ(dynG.numberOfNodes(), G.numberOfNodes());
    EXPECT_EQ(dynG.numberOfEdges(), G.numberOfEdges());
}

TEST_P(CSRGraphGTest, CtorEdgeDirection) {
    CSRGraph<int> G(10, isDirected());
    EXPECT_THAT(G.isDirected(), isDirected());
}

TEST_P(CSRGraphGTest, AddNode) {
    CSRGraph<int> G(0, isDirected());
    ASSERT_EQ(G.numberOfNodes(), 0);

    int newNode = G.addNode();
    EXPECT_EQ(newNode, 0);
    EXPECT_EQ(G.numberOfNodes(), 1);
}

INSTANTIATE_TEST_SUITE_P(TestCSRGraph, CSRGraphGTest,
                         testing::Values(GraphProperties{false}, GraphProperties{true}));
} // namespace

} // namespace NetworKit

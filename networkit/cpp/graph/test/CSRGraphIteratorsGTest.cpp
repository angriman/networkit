#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <vector>

#include <networkit/graph/CSRGraph.hpp>

namespace NetworKit {
namespace {

using ::testing::Each;
using ::testing::ElementsAreArray;
using ::testing::UnorderedElementsAreArray;

struct GraphProperties {
    bool isDirected;
};

class CSRGraphIteratorsGTest : public ::testing::TestWithParam<GraphProperties> {
protected:
    bool isDirected() const noexcept { return GetParam().isDirected; }

    template <class GraphType>
    GraphType generateCompleteGraph(int n) const {
        GraphType G(n, isDirected());
        for (int source = 0; source < n - 1; ++source) {
            for (int target = source + 1; target < n; ++target) {
                G.addEdge(source, target);
                if (isDirected()) {
                    G.addEdge(target, source);
                }
            }
        }

        return G;
    }

    template <class GraphType, class EdgeWeightType>
    GraphType generateCompleteWeightedGraph(int n, EdgeWeightType weight) const {
        GraphType G(n, isDirected());
        for (int source = 0; source < n - 1; ++source) {
            for (int target = source + 1; target < n; ++target) {
                G.addEdge(source, target, weight);
                if (isDirected()) {
                    G.addEdge(target, source, weight);
                }
            }
        }

        return G;
    }
};

TEST_P(CSRGraphIteratorsGTest, NodeIteratorIteratesOverAllNodes) {
    CSRGraph<int> G(4, isDirected());
    std::vector<int> nodes(G.nodeRange().begin(), G.nodeRange().end());
    EXPECT_THAT(nodes, ElementsAreArray({0, 1, 2, 3}));
}

TEST_P(CSRGraphIteratorsGTest, WeightedNodeIteratorIteratesOverAllNodes) {
    CSRGraph<int, float> G(4, isDirected());
    std::vector<int> nodes(G.nodeRange().begin(), G.nodeRange().end());
    EXPECT_THAT(nodes, ElementsAreArray({0, 1, 2, 3}));
}

TEST_P(CSRGraphIteratorsGTest, NeighborIteratorIteratesOverAllNeighbors) {
    auto G = generateCompleteGraph<CSRGraph<int>>(4);
    std::vector<int> neighbors(G.neighborRange(0).begin(), G.neighborRange(0).end());
    EXPECT_THAT(neighbors, UnorderedElementsAreArray({1, 2, 3}));
}

TEST_P(CSRGraphIteratorsGTest, WeightedNeighborIteratorIteratesOverAllNeighbors) {
    auto G = generateCompleteGraph<CSRGraph<int, float>>(4);
    std::vector<int> neighbors(G.neighborRange(0).begin(), G.neighborRange(0).end());
    EXPECT_THAT(neighbors, UnorderedElementsAreArray({1, 2, 3}));
}

TEST_P(CSRGraphIteratorsGTest, NeighborWeightIteratorIteratesOverAllNeigbors) {
    auto G = generateCompleteWeightedGraph<CSRGraph<int>>(4, 4);
    std::vector<int> neighbors;
    std::vector<edgeweight> weights;
    for (auto [neighbor, weight] : G.neighborWeightRange(0)) {
        neighbors.push_back(neighbor);
        weights.push_back(weight);
    }
    EXPECT_THAT(neighbors, UnorderedElementsAreArray({1, 2, 3}));
    EXPECT_THAT(weights, UnorderedElementsAreArray(
                             {defaultEdgeWeight, defaultEdgeWeight, defaultEdgeWeight}));
}

TEST_P(CSRGraphIteratorsGTest, WeightedNeighborWeightIteratorIteratesOverAllNeigbors) {
    auto G = generateCompleteWeightedGraph<CSRGraph<int, float>>(4, 4.f);
    std::vector<int> neighbors;
    std::vector<float> weights;
    for (auto [neighbor, weight] : G.neighborWeightRange(0)) {
        neighbors.push_back(neighbor);
        weights.push_back(weight);
    }
    EXPECT_THAT(neighbors, UnorderedElementsAreArray({1, 2, 3}));
    EXPECT_THAT(weights, UnorderedElementsAreArray({4.f, 4.f, 4.f}));
}

TEST_P(CSRGraphIteratorsGTest, NodeIteratorIncrement) {
    CSRGraph<int> G(4, isDirected());
    auto it = G.nodeRange().begin();

    ASSERT_EQ(*it, 0);

    EXPECT_EQ(*(++it), 1);
    EXPECT_EQ(*(it++), 1);
    EXPECT_EQ(*it, 2);

    EXPECT_EQ(*(--it), 1);
    EXPECT_EQ(*(it--), 1);
    EXPECT_EQ(*it, 0);
}

TEST_P(CSRGraphIteratorsGTest, WeightedNodeIteratorIncrement) {
    CSRGraph<int, float> G(4, isDirected());
    auto it = G.nodeRange().begin();

    ASSERT_EQ(*it, 0);

    EXPECT_EQ(*(++it), 1);
    EXPECT_EQ(*(it++), 1);
    EXPECT_EQ(*it, 2);

    EXPECT_EQ(*(--it), 1);
    EXPECT_EQ(*(it--), 1);
    EXPECT_EQ(*it, 0);
}

TEST_P(CSRGraphIteratorsGTest, ParallelForNodesCoversAllNodesOnce) {
    CSRGraph<int> G(10, isDirected());

    std::vector<int> nodeCounts(G.upperNodeIdBound(), 0);

    G.parallelForNodes([&nodeCounts](int u) {
#pragma omp atomic
        ++nodeCounts[u];
    });

    EXPECT_THAT(nodeCounts, Each(1));
}

TEST_P(CSRGraphIteratorsGTest, ParallelForNodesYieldsCorrectNodes) {
    CSRGraph<int> G(8, isDirected());
    std::vector<int> yieldedNodes(8);

    G.parallelForNodes([&yieldedNodes](int u) { yieldedNodes[u] = u; });

    EXPECT_THAT(yieldedNodes, UnorderedElementsAreArray({0, 1, 2, 3, 4, 5, 6, 7}));
}

INSTANTIATE_TEST_SUITE_P(TestCSRGraphIterators, CSRGraphIteratorsGTest,
                         testing::Values(GraphProperties{false}, GraphProperties{true}));

} // namespace
} // namespace NetworKit

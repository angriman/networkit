#include <algorithm>
#include <atomic>
#include <concepts>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/CSRGraph.hpp>
#include <networkit/graph/test/GraphTestUtils.hpp>

namespace NetworKit {
namespace {

using ::testing::ElementsAre;
using ::testing::Pair;

template <class TestT>
class CSRGraphGTest : public ::testing::Test {
public:
    using NodeT = typename TestT::NodeT;
    using EdgeWeightT = typename TestT::EdgeWeightT;
};

TYPED_TEST_SUITE_P(CSRGraphGTest);

TYPED_TEST_P(CSRGraphGTest, testDefaultConstructor) {
    const CSRGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT> G{};
    EXPECT_EQ(G.numberOfNodes(), 0);
    EXPECT_EQ(G.upperNodeIdBound(), 0);
    EXPECT_EQ(G.numberOfEdges(), 0);
    EXPECT_NE(G.isWeighted(), (std::is_same_v<typename TestFixture::EdgeWeightT, Unweighted>));
    EXPECT_FALSE(G.isDirected());
}

TYPED_TEST_P(CSRGraphGTest, testIsDirectedConstructor) {
    const CSRGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT> G(
        TypeParam::directed);
    EXPECT_EQ(G.numberOfNodes(), 0);
    EXPECT_EQ(G.upperNodeIdBound(), 0);
    EXPECT_EQ(G.numberOfEdges(), 0);
    EXPECT_NE(G.isWeighted(), (std::is_same_v<typename TestFixture::EdgeWeightT, Unweighted>));
    EXPECT_EQ(G.isDirected(), TypeParam::directed);
}

TYPED_TEST_P(CSRGraphGTest, testNumNodesConstructor) {
    constexpr size_t kNumNodes = 10;
    const CSRGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT> G(
        kNumNodes, TypeParam::directed);
    EXPECT_EQ(G.numberOfNodes(), kNumNodes);
    EXPECT_EQ(G.upperNodeIdBound(), kNumNodes);
    EXPECT_EQ(G.numberOfEdges(), 0);
    EXPECT_NE(G.isWeighted(), (std::is_same_v<typename TestFixture::EdgeWeightT, Unweighted>));
    EXPECT_EQ(G.isDirected(), TypeParam::directed);
}

TYPED_TEST_P(CSRGraphGTest, testHasNode) {
    constexpr size_t kNumNodes = 2;
    using NodeT = typename TestFixture::NodeT;
    CSRGraph<NodeT, typename TestFixture::EdgeWeightT> G(kNumNodes, TypeParam::directed);
    EXPECT_TRUE(G.hasNode(NodeT{0}));
    EXPECT_TRUE(G.hasNode(NodeT{1}));
    EXPECT_FALSE(G.hasNode(NodeT{2}));
    EXPECT_FALSE(G.hasNode(NodeT{5}));
}

TYPED_TEST_P(CSRGraphGTest, testAddEdge) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    constexpr size_t kNumNodes = 4;
    CSRGraph<NodeT, EdgeWeightT> G(kNumNodes, TypeParam::directed);

    // addEdge reports success. A weight can always be supplied; it is accepted but ignored for
    // unweighted graphs.
    EXPECT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    EXPECT_TRUE(G.addEdge(NodeT{1}, NodeT{2}));
    EXPECT_EQ(G.numberOfEdges(), 2);
    EXPECT_EQ(G.numberOfSelfLoops(), 0);

    // A self-loop is stored only once and counted once, regardless of directedness.
    EXPECT_TRUE(G.addEdge(NodeT{3}, NodeT{3}, 5.0));
    EXPECT_EQ(G.numberOfEdges(), 3);
    EXPECT_EQ(G.numberOfSelfLoops(), 1);

    // Parallel edges are permitted: re-adding an existing edge succeeds and grows the edge count.
    EXPECT_TRUE(G.addEdge(NodeT{0}, NodeT{1}));
    EXPECT_EQ(G.numberOfEdges(), 4);
}

TYPED_TEST_P(CSRGraphGTest, testNodeRange) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    constexpr size_t kNumNodes = 4;
    const GraphT G(kNumNodes, TypeParam::directed);

    // Nodes are contiguous in [0, numberOfNodes()), so the range yields every node in order.
    const std::vector<NodeT> nodes(G.nodeRange().begin(), G.nodeRange().end());
    EXPECT_THAT(nodes, ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{3}));

    // An empty graph yields an empty node range.
    const GraphT empty{};
    EXPECT_EQ(empty.nodeRange().begin(), empty.nodeRange().end());
}

TYPED_TEST_P(CSRGraphGTest, testEdgeRange) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    CSRGraph<NodeT, EdgeWeightT> G(4, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, 2.0);
    G.addEdge(NodeT{1}, NodeT{2}, 3.0);

    // Undirected edges are reported once (with u <= v); directed edges follow the stored direction.
    // Both orderings coincide for this graph.
    const std::vector<EdgeT<NodeT>> edges(G.edgeRange().begin(), G.edgeRange().end());
    EXPECT_THAT(edges, ElementsAre(EdgeEq(NodeT{0}, NodeT{1}), EdgeEq(NodeT{1}, NodeT{2})));
}

TYPED_TEST_P(CSRGraphGTest, testEdgeWeightRange) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    GraphT G(4, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, 2.0);
    G.addEdge(NodeT{1}, NodeT{2}, 3.0);

    const std::vector<WeightedEdgeT<NodeT, EdgeWeightT>> edges(G.edgeWeightRange().begin(),
                                                               G.edgeWeightRange().end());
    ASSERT_EQ(edges.size(), 2u);
    EXPECT_THAT(edges, ElementsAre(EdgeEq(NodeT{0}, NodeT{1}), EdgeEq(NodeT{1}, NodeT{2})));

    // Weights are only meaningful for weighted graphs; unweighted graphs carry an empty weight.
    if constexpr (GraphT::isWeighted()) {
        EXPECT_EQ(edges[0].weight, EdgeWeightT{2});
        EXPECT_EQ(edges[1].weight, EdgeWeightT{3});
    }
}

TYPED_TEST_P(CSRGraphGTest, testNeighborRange) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    CSRGraph<NodeT, EdgeWeightT> G(4, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, 2.0);
    G.addEdge(NodeT{0}, NodeT{2}, 3.0);

    const std::vector<NodeT> neighbors(G.neighborRange(NodeT{0}).begin(),
                                       G.neighborRange(NodeT{0}).end());
    EXPECT_THAT(neighbors, ElementsAre(NodeT{1}, NodeT{2}));

    // For undirected graphs the reverse edge is stored too, so node 1 (whose adjacency segment is
    // not the first one) lists node 0 as its neighbor. Directed graphs only store the forward edge.
    const std::vector<NodeT> neighborsOfOne(G.neighborRange(NodeT{1}).begin(),
                                            G.neighborRange(NodeT{1}).end());
    if (TypeParam::directed) {
        EXPECT_TRUE(neighborsOfOne.empty());
    } else {
        EXPECT_THAT(neighborsOfOne, ElementsAre(NodeT{0}));
    }

    // A node without outgoing edges yields an empty range.
    EXPECT_EQ(G.neighborRange(NodeT{3}).begin(), G.neighborRange(NodeT{3}).end());
}

TYPED_TEST_P(CSRGraphGTest, testNeighborIteratorSemantics) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    GraphT G(3, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, 2.0);
    G.addEdge(NodeT{0}, NodeT{2}, 3.0);

    const auto range = G.neighborRange(NodeT{0});
    auto it = range.begin();
    EXPECT_EQ(*it, NodeT{1});

    // Post-increment yields the pre-increment position, then advances.
    const auto prev = it++;
    EXPECT_EQ(*prev, NodeT{1});
    EXPECT_EQ(*it, NodeT{2});

    // Pre-increment reaches the end and compares equal to end().
    EXPECT_NE(it, range.end());
    EXPECT_EQ(++it, range.end());

    // Iterator is a forward iterator: copies advance independently of the original.
    auto a = range.begin();
    auto b = a;
    ++a;
    EXPECT_EQ(*b, NodeT{1});
    EXPECT_NE(a, b);
    ++b;
    EXPECT_EQ(a, b);

    // The default-constructed iterator (required by the Python bindings) is constructible.
    typename GraphT::NeighborIterator def;
    (void)def;
}

TYPED_TEST_P(CSRGraphGTest, testNeighborWeightIteratorSemantics) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    GraphT G(3, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, 2.0);
    G.addEdge(NodeT{0}, NodeT{2}, 3.0);

    const auto range = G.weightNeighborRange(NodeT{0});
    auto it = range.begin();
    EXPECT_EQ((*it).first, NodeT{1});

    // Post-increment yields the pre-increment position, then advances.
    const auto prev = it++;
    EXPECT_EQ((*prev).first, NodeT{1});
    EXPECT_EQ((*it).first, NodeT{2});

    // Pre-increment reaches the end and compares equal to end().
    EXPECT_NE(it, range.end());
    EXPECT_EQ(++it, range.end());

    // Iterator is a forward iterator: copies advance independently of the original.
    auto a = range.begin();
    auto b = a;
    ++a;
    EXPECT_EQ((*b).first, NodeT{1});
    EXPECT_NE(a, b);
    ++b;
    EXPECT_EQ(a, b);

    // The default-constructed iterator (required by the Python bindings) is constructible.
    typename GraphT::NeighborWeightIterator def;
    (void)def;
}

TYPED_TEST_P(CSRGraphGTest, testNeighborWeightRange) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    GraphT G(4, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, 2.0);
    G.addEdge(NodeT{0}, NodeT{2}, 3.0);

    const std::vector<std::pair<NodeT, EdgeWeightT>> neighbors(
        G.weightNeighborRange(NodeT{0}).begin(), G.weightNeighborRange(NodeT{0}).end());
    ASSERT_EQ(neighbors.size(), 2u);
    EXPECT_EQ(neighbors[0].first, NodeT{1});
    EXPECT_EQ(neighbors[1].first, NodeT{2});

    // Weights are only meaningful for weighted graphs; unweighted graphs carry an empty weight.
    if constexpr (GraphT::isWeighted()) {
        EXPECT_EQ(neighbors[0].second, EdgeWeightT{2});
        EXPECT_EQ(neighbors[1].second, EdgeWeightT{3});
    }
}

TYPED_TEST_P(CSRGraphGTest, testForNodes) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const CSRGraph<NodeT, EdgeWeightT> G(4, TypeParam::directed);

    std::vector<NodeT> visited;
    G.forNodes([&](NodeT u) { visited.push_back(u); });
    EXPECT_THAT(visited, ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{3}));
}

TYPED_TEST_P(CSRGraphGTest, testParallelForNodes) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    constexpr size_t kNumNodes = 100;
    const CSRGraph<NodeT, EdgeWeightT> G(kNumNodes, TypeParam::directed);

    // Each node is visited exactly once; writing to a distinct slot per node avoids data races.
    std::vector<char> visited(kNumNodes, 0);
    G.parallelForNodes([&](NodeT u) { visited[u] = 1; });
    EXPECT_EQ(static_cast<size_t>(std::count(visited.begin(), visited.end(), 1)), kNumNodes);
}

TYPED_TEST_P(CSRGraphGTest, testForEdges) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    GraphT G(4, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, 2.0);
    G.addEdge(NodeT{0}, NodeT{2}, 3.0);

    // The (u, v) overload. Undirected edges are reported once with u <= v, so the reported set
    // coincides with the directed case for this graph.
    std::vector<EdgeT<NodeT>> edges;
    G.forEdges([&](NodeT u, NodeT v) { edges.emplace_back(u, v); });
    EXPECT_THAT(edges, ElementsAre(EdgeEq(NodeT{0}, NodeT{1}), EdgeEq(NodeT{0}, NodeT{2})));

    // The (u, v, weight) overload is only meaningful for weighted graphs.
    if constexpr (GraphT::isWeighted()) {
        std::vector<WeightedEdgeT<NodeT, EdgeWeightT>> weighted;
        G.forEdges([&](NodeT u, NodeT v, EdgeWeightT w) { weighted.emplace_back(u, v, w); });
        EXPECT_THAT(weighted, ElementsAre(WeightedEdgeEq(NodeT{0}, NodeT{1}, EdgeWeightT{2}),
                                          WeightedEdgeEq(NodeT{0}, NodeT{2}, EdgeWeightT{3})));
    }
}

TYPED_TEST_P(CSRGraphGTest, testParallelForEdges) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    CSRGraph<NodeT, EdgeWeightT> G(4, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1});
    G.addEdge(NodeT{0}, NodeT{2});
    G.addEdge(NodeT{1}, NodeT{2});

    // Every edge is reported exactly once, so the count matches numberOfEdges().
    std::atomic<count> edgeCount{0};
    G.parallelForEdges([&](NodeT, NodeT) { ++edgeCount; });
    EXPECT_EQ(edgeCount.load(), G.numberOfEdges());
}

TYPED_TEST_P(CSRGraphGTest, testForNeighborsOf) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    GraphT G(4, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, 2.0);
    G.addEdge(NodeT{0}, NodeT{2}, 3.0);

    std::vector<NodeT> neighbors;
    G.forNeighborsOf(NodeT{0}, [&](NodeT v) { neighbors.push_back(v); });
    EXPECT_THAT(neighbors, ElementsAre(NodeT{1}, NodeT{2}));

    // A node without outgoing edges yields no callbacks.
    G.forNeighborsOf(NodeT{3},
                     [&](NodeT) { ADD_FAILURE() << "unexpected neighbor of isolated node"; });

    // The (v, weight) overload is only meaningful for weighted graphs.
    if constexpr (GraphT::isWeighted()) {
        std::vector<std::pair<NodeT, EdgeWeightT>> weighted;
        G.forNeighborsOf(NodeT{0}, [&](NodeT v, EdgeWeightT w) { weighted.emplace_back(v, w); });
        EXPECT_THAT(weighted,
                    ElementsAre(Pair(NodeT{1}, EdgeWeightT{2}), Pair(NodeT{2}, EdgeWeightT{3})));
    }
}

TYPED_TEST_P(CSRGraphGTest, testForEdgesOf) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    CSRGraph<NodeT, EdgeWeightT> G(4, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, 2.0);
    G.addEdge(NodeT{0}, NodeT{2}, 3.0);

    // forEdgesOf reports every incident (outgoing) edge; the source is always the queried node.
    std::vector<EdgeT<NodeT>> edges;
    G.forEdgesOf(NodeT{0}, [&](NodeT u, NodeT v) { edges.emplace_back(u, v); });
    EXPECT_THAT(edges, ElementsAre(EdgeEq(NodeT{0}, NodeT{1}), EdgeEq(NodeT{0}, NodeT{2})));
}

REGISTER_TYPED_TEST_SUITE_P(CSRGraphGTest, testDefaultConstructor, testIsDirectedConstructor,
                            testNumNodesConstructor, testHasNode, testAddEdge, testNodeRange,
                            testEdgeRange, testEdgeWeightRange, testNeighborRange,
                            testNeighborWeightRange, testNeighborIteratorSemantics,
                            testNeighborWeightIteratorSemantics, testForNodes, testParallelForNodes,
                            testForEdges, testParallelForEdges, testForNeighborsOf, testForEdgesOf);

template <class NodeType, class EdgeWeightType, bool Directed>
struct CSRGraphConfig {
    using NodeT = NodeType;
    using EdgeWeightT = EdgeWeightType;
    static constexpr bool directed = Directed;
};

using CSRGraphTestTypes =
    ::testing::Types<CSRGraphConfig<node, edgeweight, false>,
                     CSRGraphConfig<node, edgeweight, true>, CSRGraphConfig<int, float, false>,
                     CSRGraphConfig<int, float, true>, CSRGraphConfig<int, Unweighted, false>,
                     CSRGraphConfig<int, Unweighted, true>>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestCSRGraph, CSRGraphGTest, CSRGraphTestTypes, );

} // namespace
} // namespace NetworKit

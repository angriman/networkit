#include <algorithm>
#include <atomic>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <forward_list>
#include <initializer_list>
#include <random>
#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/CSRGraph.hpp>
#include <networkit/graph/test/GraphTestUtils.hpp>

namespace NetworKit {
namespace {

using ::testing::_;
using ::testing::ElementsAre;
using ::testing::Pair;

template <class TestT>
class CSRGraphGTest : public ::testing::Test {
public:
    using NodeT = typename TestT::NodeT;
    using EdgeWeightT = typename TestT::EdgeWeightT;
};

// Records the arguments that each supported handle signature of the neighborhood and edge
// iterators receives, with the weight (0 if the handle takes none) converted to double.
template <class GraphT>
std::vector<std::tuple<typename GraphT::NodeT, typename GraphT::NodeT, double>>
recordHandleArguments(const GraphT &G) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
    // AdjListGraph::isWeighted() is not constexpr, so look at the weight type instead.
    constexpr bool kWeighted = !std::same_as<EdgeWeightT, Unweighted>;
    std::vector<std::tuple<NodeT, NodeT, double>> args;
    const auto recordEdge = [&args](NodeT x, NodeT y) -> void { args.emplace_back(x, y, 0.); };
    const auto recordNeighbor = [&args](NodeT v) -> void { args.emplace_back(v, v, 0.); };

    G.forNodes([&](NodeT u) -> void {
        G.forNeighborsOf(u, recordNeighbor);
        G.forNeighborsOf(u, recordEdge);
        G.forEdgesOf(u, recordNeighbor);
        G.forEdgesOf(u, recordEdge);
        if constexpr (kWeighted) {
            G.forNeighborsOf(
                u, [&args](NodeT v, EdgeWeightT w) -> void { args.emplace_back(v, v, w); });
            G.forNeighborsOf(u, [&args](NodeT x, NodeT y, EdgeWeightT w) -> void {
                args.emplace_back(x, y, w);
            });
            G.forEdgesOf(u, [&args](NodeT x, NodeT y, EdgeWeightT w) -> void {
                args.emplace_back(x, y, w);
            });
        }
    });
    G.forEdges(recordEdge);
    if constexpr (kWeighted) {
        G.forEdges(
            [&args](NodeT x, NodeT y, EdgeWeightT w) -> void { args.emplace_back(x, y, w); });
    }
    return args;
}

// The (neighbor, weight) pairs of every node of `G` in iteration order, with the weight converted
// to double (0 for unweighted graphs).
template <class GraphT>
std::vector<std::vector<std::pair<typename GraphT::NodeT, double>>>
weightedNeighborhoods(const GraphT &G) {
    using NodeT = typename GraphT::NodeT;
    std::vector<std::vector<std::pair<NodeT, double>>> neighborhoods(G.upperNodeIdBound());
    G.forNodes([&](NodeT u) -> void {
        for (const auto &[v, w] : G.weightNeighborRange(u)) {
            if constexpr (GraphT::isWeighted()) {
                neighborhoods[u].emplace_back(v, w);
            } else {
                neighborhoods[u].emplace_back(v, 0.);
            }
        }
    });
    return neighborhoods;
}

TYPED_TEST_SUITE_P(CSRGraphGTest);

TYPED_TEST_P(CSRGraphGTest, testDefaultConstructor) {
    const CSRGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT> G{};
    EXPECT_EQ(G.numberOfNodes(), 0);
    EXPECT_EQ(G.upperNodeIdBound(), 0);
    EXPECT_EQ(G.numberOfEdges(), 0);
    EXPECT_NE(G.isWeighted(), (std::is_same_v<typename TestFixture::EdgeWeightT, Unweighted>));
    EXPECT_FALSE(G.isDirected());
}

TYPED_TEST_P(CSRGraphGTest, testTypeAliases) {
    using GraphT = CSRGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT>;
    // Algorithms templated on the graph type (e.g. BFS) read the node and weight types from these.
    static_assert(std::same_as<typename GraphT::NodeT, typename TestFixture::NodeT>);
    static_assert(std::same_as<typename GraphT::EdgeWeightT, typename TestFixture::EdgeWeightT>);
}

TYPED_TEST_P(CSRGraphGTest, testNeighborIteratorConcepts) {
    using GraphT = CSRGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT>;
    static_assert(std::forward_iterator<typename GraphT::NeighborIterator>);
    static_assert(std::forward_iterator<typename GraphT::NeighborWeightIterator>);
    static_assert(std::ranges::forward_range<typename GraphT::NeighborRange>);
    static_assert(std::ranges::forward_range<typename GraphT::NeighborWeightRange>);
    // Default-constructed iterators compare equal.
    EXPECT_EQ(typename GraphT::NeighborIterator{}, typename GraphT::NeighborIterator{});
    EXPECT_EQ(typename GraphT::NeighborWeightIterator{}, typename GraphT::NeighborWeightIterator{});
}

TYPED_TEST_P(CSRGraphGTest, testSingleArgumentConstructor) {
    using NodeT = typename TestFixture::NodeT;
    using GraphT = CSRGraph<NodeT, typename TestFixture::EdgeWeightT>;
    // The number of nodes can be passed as any integer type without ambiguity. The graph is
    // undirected by default.
    const GraphT fromLiteral(4);
    const GraphT fromCount(count{4});
    const GraphT fromSizeT(std::size_t{4});
    const GraphT fromNode(NodeT{4});
    for (const GraphT *G : {&fromLiteral, &fromCount, &fromSizeT, &fromNode}) {
        EXPECT_EQ(G->numberOfNodes(), 4);
        EXPECT_EQ(G->numberOfEdges(), 0);
        EXPECT_FALSE(G->isDirected());
    }
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
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{1}, NodeT{2}));
    EXPECT_EQ(G.numberOfEdges(), 2);
    EXPECT_EQ(G.numberOfSelfLoops(), 0);

    // A self-loop is stored only once and counted once, regardless of directedness.
    ASSERT_TRUE(G.addEdge(NodeT{3}, NodeT{3}, 5.0));
    EXPECT_EQ(G.numberOfEdges(), 3);
    EXPECT_EQ(G.numberOfSelfLoops(), 1);

    // Parallel edges are permitted: re-adding an existing edge succeeds and grows the edge count.
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}));
    EXPECT_EQ(G.numberOfEdges(), 4);
}

TYPED_TEST_P(CSRGraphGTest, testAddEdgeKeepsAdjacencySorted) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    GraphT G(5, TypeParam::directed);
    // Unsorted targets, with the parallel edge (0, 3) inserted around other edges. The weight
    // identifies each edge.
    for (const auto &[u, v, w] : {std::tuple{0, 3, 1}, std::tuple{0, 1, 2}, std::tuple{0, 3, 3},
                                  std::tuple{0, 2, 4}, std::tuple{4, 3, 5}, std::tuple{1, 3, 6}}) {
        ASSERT_TRUE(G.addEdge(NodeT(u), NodeT(v), w));
    }

    // Weights move with their edges, and parallel edges stay in insertion order.
    const auto neighborhoods = weightedNeighborhoods(G);
    const double w = GraphT::isWeighted() ? 1. : 0.;
    EXPECT_THAT(neighborhoods[0], ElementsAre(Pair(NodeT{1}, 2 * w), Pair(NodeT{2}, 4 * w),
                                              Pair(NodeT{3}, 1 * w), Pair(NodeT{3}, 3 * w)));
    if (TypeParam::directed) {
        EXPECT_THAT(neighborhoods[3], ElementsAre());
    } else {
        // The reverse entries are kept sorted too.
        EXPECT_THAT(neighborhoods[3], ElementsAre(Pair(NodeT{0}, 1 * w), Pair(NodeT{0}, 3 * w),
                                                  Pair(NodeT{1}, 6 * w), Pair(NodeT{4}, 5 * w)));
    }
}

TYPED_TEST_P(CSRGraphGTest, testFromEdges) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    // Unsorted input with a self-loop and parallel edges.
    const std::vector<WeightedEdgeT<NodeT, double>> input{
        {0, 3, 4.}, {0, 1, 2.}, {2, 2, 5.}, {2, 0, 3.}, {0, 1, 6.}};
    const GraphT G = GraphT::fromEdges(4, input, TypeParam::directed);

    EXPECT_EQ(G.numberOfNodes(), 4);
    EXPECT_EQ(G.numberOfEdges(), input.size());
    EXPECT_EQ(G.numberOfSelfLoops(), 1);
    EXPECT_EQ(G.isDirected(), TypeParam::directed);

    const auto neighborsOf = [&G](NodeT u) -> std::vector<NodeT> {
        return {G.neighborRange(u).begin(), G.neighborRange(u).end()};
    };
    if (TypeParam::directed) {
        EXPECT_THAT(neighborsOf(0), ElementsAre(NodeT{1}, NodeT{1}, NodeT{3}));
        EXPECT_THAT(neighborsOf(1), ElementsAre());
        EXPECT_THAT(neighborsOf(2), ElementsAre(NodeT{0}, NodeT{2}));
        EXPECT_THAT(neighborsOf(3), ElementsAre());
    } else {
        EXPECT_THAT(neighborsOf(0), ElementsAre(NodeT{1}, NodeT{1}, NodeT{2}, NodeT{3}));
        EXPECT_THAT(neighborsOf(1), ElementsAre(NodeT{0}, NodeT{0}));
        EXPECT_THAT(neighborsOf(2), ElementsAre(NodeT{0}, NodeT{2}));
        EXPECT_THAT(neighborsOf(3), ElementsAre(NodeT{0}));
    }

    // Weights travel with their edges, and the parallel edges (0, 1) keep their input order.
    if constexpr (GraphT::isWeighted()) {
        const auto neighborhoods = weightedNeighborhoods(G);
        EXPECT_THAT(neighborhoods[2], ElementsAre(Pair(NodeT{0}, 3.), Pair(NodeT{2}, 5.)));
        if (TypeParam::directed) {
            EXPECT_THAT(neighborhoods[0],
                        ElementsAre(Pair(NodeT{1}, 2.), Pair(NodeT{1}, 6.), Pair(NodeT{3}, 4.)));
        } else {
            EXPECT_THAT(neighborhoods[0], ElementsAre(Pair(NodeT{1}, 2.), Pair(NodeT{1}, 6.),
                                                      Pair(NodeT{2}, 3.), Pair(NodeT{3}, 4.)));
            EXPECT_THAT(neighborhoods[1], ElementsAre(Pair(NodeT{0}, 2.), Pair(NodeT{0}, 6.)));
        }
    }

    // Plain edges suffice for unweighted graphs.
    if constexpr (!GraphT::isWeighted()) {
        const std::vector<EdgeT<NodeT>> plain{{0, 1}, {1, 2}};
        const GraphT H = GraphT::fromEdges(3, plain, TypeParam::directed);
        EXPECT_EQ(H.numberOfEdges(), 2);
    }
}

TYPED_TEST_P(CSRGraphGTest, testFromEdgesWithoutEdges) {
    using GraphT = CSRGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT>;
    const std::vector<WeightedEdgeT<typename TestFixture::NodeT, double>> noEdges;

    const GraphT empty = GraphT::fromEdges(0, noEdges, TypeParam::directed);
    EXPECT_EQ(empty.numberOfNodes(), 0);
    EXPECT_EQ(empty.numberOfEdges(), 0);
    EXPECT_EQ(empty.isDirected(), TypeParam::directed);

    const GraphT isolated = GraphT::fromEdges(3, noEdges, TypeParam::directed);
    EXPECT_EQ(isolated.numberOfNodes(), 3);
    EXPECT_EQ(isolated.numberOfEdges(), 0);
    isolated.forNodes(
        [&](typename TestFixture::NodeT u) -> void { EXPECT_EQ(isolated.degree(u), 0); });
}

TYPED_TEST_P(CSRGraphGTest, testFromGraphMatchesAdjListGraph) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    constexpr count kNumNodes = 5;
    AdjListGraph<NodeT, EdgeWeightT> adj(kNumNodes, GraphT::isWeighted(), TypeParam::directed);
    for (const auto &[u, v, w] : {std::tuple{3, 1, 2}, std::tuple{0, 2, 3}, std::tuple{2, 1, 4},
                                  std::tuple{4, 4, 5}, std::tuple{0, 1, 6}, std::tuple{1, 0, 7}}) {
        if constexpr (GraphT::isWeighted()) {
            adj.addEdge(NodeT(u), NodeT(v), EdgeWeightT(w));
        } else {
            adj.addEdge(NodeT(u), NodeT(v));
        }
    }

    const GraphT csr = GraphT::fromGraph(adj);
    EXPECT_EQ(csr.numberOfNodes(), adj.numberOfNodes());
    EXPECT_EQ(csr.numberOfEdges(), adj.numberOfEdges());
    EXPECT_EQ(csr.numberOfSelfLoops(), adj.numberOfSelfLoops());
    EXPECT_EQ(csr.isDirected(), adj.isDirected());

    // Same neighborhoods with the same weights; CSRGraph additionally sorts them by target, keeping
    // parallel edges in insertion order.
    adj.forNodes([&](NodeT u) -> void {
        EXPECT_EQ(csr.degree(u), adj.degree(u));
        std::vector<std::pair<NodeT, double>> expected, actual;
        if constexpr (GraphT::isWeighted()) {
            adj.forNeighborsOf(
                u, [&](NodeT v, EdgeWeightT w) -> void { expected.emplace_back(v, w); });
        } else {
            adj.forNeighborsOf(u, [&](NodeT v) -> void { expected.emplace_back(v, 0.); });
        }
        for (const auto &[v, w] : csr.weightNeighborRange(u)) {
            if constexpr (GraphT::isWeighted()) {
                actual.emplace_back(v, w);
            } else {
                actual.emplace_back(v, 0.);
            }
        }
        std::ranges::stable_sort(expected, {}, &std::pair<NodeT, double>::first);
        EXPECT_EQ(actual, expected);
    });
}

TYPED_TEST_P(CSRGraphGTest, testBuildersAgreeOnRandomGraph) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    // Few nodes and many edges, so that there are many parallel edges and self-loops.
    constexpr count kNumNodes = 50;
    constexpr count kNumEdges = 1000;
    std::mt19937 prng(42);
    std::uniform_int_distribution<int> randomNode(0, kNumNodes - 1);

    // The weight is the edge's index, which identifies it.
    std::vector<WeightedEdgeT<NodeT, double>> edgeList;
    GraphT incremental(kNumNodes, TypeParam::directed);
    AdjListGraph<NodeT, EdgeWeightT> adj(kNumNodes, GraphT::isWeighted(), TypeParam::directed);
    for (count i = 0; i < kNumEdges; ++i) {
        const auto u = NodeT(randomNode(prng));
        const auto v = NodeT(randomNode(prng));
        edgeList.emplace_back(u, v, double(i));
        ASSERT_TRUE(incremental.addEdge(u, v, double(i)));
        if constexpr (GraphT::isWeighted()) {
            adj.addEdge(u, v, EdgeWeightT(i));
        } else {
            adj.addEdge(u, v);
        }
    }

    const GraphT fromEdges = GraphT::fromEdges(kNumNodes, edgeList, TypeParam::directed);
    const GraphT fromGraph = GraphT::fromGraph(adj);
    EXPECT_EQ(fromEdges.numberOfEdges(), kNumEdges);
    EXPECT_EQ(fromEdges.numberOfSelfLoops(), adj.numberOfSelfLoops());
    for (const GraphT *G : std::initializer_list<const GraphT *>{&incremental, &fromGraph}) {
        EXPECT_EQ(G->numberOfEdges(), fromEdges.numberOfEdges());
        EXPECT_EQ(G->numberOfSelfLoops(), fromEdges.numberOfSelfLoops());
    }

    // All builders keep neighborhoods sorted by target and parallel edges in insertion order.
    const auto neighborhoods = weightedNeighborhoods(fromEdges);
    EXPECT_EQ(weightedNeighborhoods(incremental), neighborhoods);
    EXPECT_EQ(weightedNeighborhoods(fromGraph), neighborhoods);
    for (const auto &neighborhood : neighborhoods) {
        EXPECT_TRUE(std::ranges::is_sorted(neighborhood));
    }
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
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{1}, NodeT{2}, 3.0));

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
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{1}, NodeT{2}, 3.0));

    const std::vector<WeightedEdgeT<NodeT, EdgeWeightT>> edges(G.edgeWeightRange().begin(),
                                                               G.edgeWeightRange().end());
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
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{2}, 3.0));

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
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{2}, 3.0));

    const auto range = G.neighborRange(NodeT{0});
    auto it = range.begin();
    ASSERT_EQ(*it, NodeT{1});

    // Post-increment yields the pre-increment position, then advances.
    const auto prev = it++;
    ASSERT_EQ(*prev, NodeT{1});
    ASSERT_EQ(*it, NodeT{2});

    // Pre-increment reaches the end and compares equal to end().
    EXPECT_NE(it, range.end());
    EXPECT_EQ(++it, range.end());

    // Iterator is a forward iterator: copies advance independently of the original.
    auto a = range.begin();
    auto b = a;
    ++a;
    ASSERT_EQ(*b, NodeT{1});
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
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{2}, 3.0));

    const auto range = G.weightNeighborRange(NodeT{0});
    auto it = range.begin();
    ASSERT_EQ((*it).first, NodeT{1});

    // Post-increment yields the pre-increment position, then advances.
    const auto prev = it++;
    ASSERT_EQ((*prev).first, NodeT{1});
    ASSERT_EQ((*it).first, NodeT{2});

    // Pre-increment reaches the end and compares equal to end().
    EXPECT_NE(it, range.end());
    EXPECT_EQ(++it, range.end());

    // Iterator is a forward iterator: copies advance independently of the original.
    auto a = range.begin();
    auto b = a;
    ++a;
    ASSERT_EQ((*b).first, NodeT{1});
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
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{2}, 3.0));

    const std::vector<std::pair<NodeT, EdgeWeightT>> neighbors(
        G.weightNeighborRange(NodeT{0}).begin(), G.weightNeighborRange(NodeT{0}).end());
    EXPECT_THAT(neighbors, ElementsAre(Pair(NodeT{1}, _), Pair(NodeT{2}, _)));

    // Weights are only meaningful for weighted graphs; unweighted graphs carry an empty weight.
    if constexpr (GraphT::isWeighted()) {
        EXPECT_THAT(neighbors, ElementsAre(Pair(_, EdgeWeightT{2}), Pair(_, EdgeWeightT{3})));
    }
}

TYPED_TEST_P(CSRGraphGTest, testForNodes) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const CSRGraph<NodeT, EdgeWeightT> G(4, TypeParam::directed);

    std::vector<NodeT> visited;
    G.forNodes([&visited](NodeT u) { visited.push_back(u); });
    EXPECT_THAT(visited, ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{3}));
}

TYPED_TEST_P(CSRGraphGTest, testParallelForNodes) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    constexpr size_t kNumNodes = 100;
    const CSRGraph<NodeT, EdgeWeightT> G(kNumNodes, TypeParam::directed);

    // Each node is visited exactly once; writing to a distinct slot per node avoids data races.
    std::vector<char> visited(kNumNodes, 0);
    G.parallelForNodes([&visited](NodeT u) { visited[u] = 1; });
    EXPECT_EQ(static_cast<size_t>(std::ranges::count(visited, 1)), kNumNodes);
}

TYPED_TEST_P(CSRGraphGTest, testForEdges) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    GraphT G(4, TypeParam::directed);
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{2}, 3.0));

    // The (u, v) overload. Directed edges follow the stored direction, undirected edges are
    // reported once with u >= v (as in AdjListGraph).
    std::vector<EdgeT<NodeT>> edges;
    G.forEdges([&edges](NodeT u, NodeT v) -> void { edges.emplace_back(u, v); });
    if (TypeParam::directed) {
        EXPECT_THAT(edges, ElementsAre(EdgeEq(NodeT{0}, NodeT{1}), EdgeEq(NodeT{0}, NodeT{2})));
    } else {
        EXPECT_THAT(edges, ElementsAre(EdgeEq(NodeT{1}, NodeT{0}), EdgeEq(NodeT{2}, NodeT{0})));
    }

    // The (u, v, weight) overload is only meaningful for weighted graphs.
    if constexpr (GraphT::isWeighted()) {
        std::vector<WeightedEdgeT<NodeT, EdgeWeightT>> weighted;
        G.forEdges([&weighted](NodeT u, NodeT v, EdgeWeightT w) -> void {
            weighted.emplace_back(u, v, w);
        });
        if (TypeParam::directed) {
            EXPECT_THAT(weighted, ElementsAre(WeightedEdgeEq(NodeT{0}, NodeT{1}, EdgeWeightT{2}),
                                              WeightedEdgeEq(NodeT{0}, NodeT{2}, EdgeWeightT{3})));
        } else {
            EXPECT_THAT(weighted, ElementsAre(WeightedEdgeEq(NodeT{1}, NodeT{0}, EdgeWeightT{2}),
                                              WeightedEdgeEq(NodeT{2}, NodeT{0}, EdgeWeightT{3})));
        }
    }
}

TYPED_TEST_P(CSRGraphGTest, testParallelForEdges) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    CSRGraph<NodeT, EdgeWeightT> G(4, TypeParam::directed);
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}));
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{2}));
    ASSERT_TRUE(G.addEdge(NodeT{1}, NodeT{2}));

    // Every edge is reported exactly once, so the count matches numberOfEdges().
    std::atomic<count> edgeCount{0};
    G.parallelForEdges([&edgeCount](NodeT, NodeT) { ++edgeCount; });
    EXPECT_EQ(edgeCount.load(), G.numberOfEdges());
}

TYPED_TEST_P(CSRGraphGTest, testForNeighborsOf) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    GraphT G(4, TypeParam::directed);
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{2}, 3.0));

    std::vector<NodeT> neighbors;
    G.forNeighborsOf(NodeT{0}, [&neighbors](NodeT v) { neighbors.push_back(v); });
    EXPECT_THAT(neighbors, ElementsAre(NodeT{1}, NodeT{2}));

    // A two-node handle receives (u, v), as in AdjListGraph, not (v, weight) converted to nodes.
    std::vector<EdgeT<NodeT>> edges;
    G.forNeighborsOf(NodeT{0}, [&edges](NodeT u, NodeT v) -> void { edges.emplace_back(u, v); });
    EXPECT_THAT(edges, ElementsAre(EdgeEq(NodeT{0}, NodeT{1}), EdgeEq(NodeT{0}, NodeT{2})));

    // A node without outgoing edges yields no callbacks.
    G.forNeighborsOf(NodeT{3},
                     [](NodeT) { ADD_FAILURE() << "unexpected neighbor of isolated node"; });

    // The (v, weight) overload is only meaningful for weighted graphs.
    if constexpr (GraphT::isWeighted()) {
        std::vector<std::pair<NodeT, EdgeWeightT>> weighted;
        G.forNeighborsOf(NodeT{0},
                         [&weighted](NodeT v, EdgeWeightT w) { weighted.emplace_back(v, w); });
        EXPECT_THAT(weighted,
                    ElementsAre(Pair(NodeT{1}, EdgeWeightT{2}), Pair(NodeT{2}, EdgeWeightT{3})));
    }
}

TYPED_TEST_P(CSRGraphGTest, testForEdgesOf) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    CSRGraph<NodeT, EdgeWeightT> G(4, TypeParam::directed);
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{1}, 2.0));
    ASSERT_TRUE(G.addEdge(NodeT{0}, NodeT{2}, 3.0));

    // forEdgesOf reports every incident (outgoing) edge; the source is always the queried node.
    std::vector<EdgeT<NodeT>> edges;
    G.forEdgesOf(NodeT{0}, [&edges](NodeT u, NodeT v) { edges.emplace_back(u, v); });
    EXPECT_THAT(edges, ElementsAre(EdgeEq(NodeT{0}, NodeT{1}), EdgeEq(NodeT{0}, NodeT{2})));

    // A single-node handle receives the neighbor, as in AdjListGraph.
    std::vector<NodeT> neighbors;
    G.forEdgesOf(NodeT{0}, [&neighbors](NodeT v) -> void { neighbors.push_back(v); });
    EXPECT_THAT(neighbors, ElementsAre(NodeT{1}, NodeT{2}));
}

TYPED_TEST_P(CSRGraphGTest, testHandleArgumentsMatchAdjListGraph) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    using GraphT = CSRGraph<NodeT, EdgeWeightT>;
    constexpr count kNumNodes = 4;
    GraphT csr(kNumNodes, TypeParam::directed);
    AdjListGraph<NodeT, EdgeWeightT> adj(kNumNodes, GraphT::isWeighted(), TypeParam::directed);

    // Includes a self-loop, and edges whose source has a larger id than their target.
    for (const auto &[u, v, w] :
         {std::tuple{0, 1, 2}, std::tuple{0, 2, 3}, std::tuple{2, 1, 4}, std::tuple{3, 3, 5}}) {
        ASSERT_TRUE(csr.addEdge(NodeT(u), NodeT(v), w));
        if constexpr (GraphT::isWeighted()) {
            adj.addEdge(NodeT(u), NodeT(v), EdgeWeightT(w));
        } else {
            adj.addEdge(NodeT(u), NodeT(v));
        }
    }

    // Every handle signature both graphs support must receive the same arguments in the same
    // order, so that CSRGraph can replace AdjListGraph in algorithms. CSRGraph sorts each
    // neighborhood by target while AdjListGraph keeps insertion order, so the edges above are
    // chosen such that every AdjListGraph neighborhood is sorted as well.
    EXPECT_EQ(recordHandleArguments(csr), recordHandleArguments(adj));
}

REGISTER_TYPED_TEST_SUITE_P(
    CSRGraphGTest, testDefaultConstructor, testTypeAliases, testNeighborIteratorConcepts,
    testSingleArgumentConstructor, testNumNodesConstructor, testHasNode, testAddEdge,
    testAddEdgeKeepsAdjacencySorted, testFromEdges, testFromEdgesWithoutEdges,
    testFromGraphMatchesAdjListGraph, testBuildersAgreeOnRandomGraph, testNodeRange, testEdgeRange,
    testEdgeWeightRange, testNeighborRange, testNeighborWeightRange, testNeighborIteratorSemantics,
    testNeighborWeightIteratorSemantics, testForNodes, testParallelForNodes, testForEdges,
    testParallelForEdges, testForNeighborsOf, testForEdgesOf, testHandleArgumentsMatchAdjListGraph);

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

TEST(CSRGraphHandleGTest, testWeightTypeSameAsEdgeId) {
    // CSRGraph has no edge ids, so a (u, v, weight) handle receives the weight even when the
    // weight type is the same as edgeid, which would otherwise select the (u, v, edgeId) signature.
    static_assert(std::same_as<uint64_t, edgeid>);
    CSRGraph<node, uint64_t> G(2);
    ASSERT_TRUE(G.addEdge(0, 1, 7));

    std::vector<std::tuple<node, node, uint64_t>> edges;
    G.forEdgesOf(0, [&edges](node u, node v, uint64_t w) -> void { edges.emplace_back(u, v, w); });
    EXPECT_THAT(edges, ElementsAre(std::tuple{node{0}, node{1}, uint64_t{7}}));
}

// Whether CSRGraph::fromEdges() accepts an edge list of the given type.
template <class GraphT, class EdgeList>
concept BuildableFromEdges =
    requires(const EdgeList &edgeList) { GraphT::fromEdges(count{0}, edgeList); };

TEST(CSRGraphFromEdgesGTest, testEdgeListConstraints) {
    using Weighted = CSRGraph<int, float>;
    using Plain = CSRGraph<int, Unweighted>;
    static_assert(BuildableFromEdges<Weighted, std::vector<WeightedEdgeT<int, double>>>);
    static_assert(BuildableFromEdges<Plain, std::vector<WeightedEdgeT<int, double>>>);
    static_assert(BuildableFromEdges<Plain, std::vector<EdgeT<int>>>);
    // Weighted graphs need weights, and all graphs need endpoints.
    static_assert(!BuildableFromEdges<Weighted, std::vector<EdgeT<int>>>);
    static_assert(!BuildableFromEdges<Plain, std::vector<int>>);
}

TEST(CSRGraphFromEdgesGTest, testForwardOnlyEdgeList) {
    const std::forward_list<WeightedEdgeT<node, double>> edgeList{{1, 0, 2.}, {0, 2, 3.}};
    const auto G = CSRGraph<node, edgeweight>::fromEdges(3, edgeList);
    EXPECT_EQ(G.numberOfEdges(), 2);
    EXPECT_THAT(weightedNeighborhoods(G)[0], ElementsAre(Pair(node{1}, 2.), Pair(node{2}, 3.)));
}

TEST(CSRGraphFromEdgesGTest, testEdgeListWithoutConstIteration) {
    // Like std::views::filter, this range can only be iterated when it is not const.
    struct MutableRange {
        std::vector<EdgeT<node>> edgeList;
        auto begin() { return edgeList.begin(); }
        auto end() { return edgeList.end(); }
    };
    MutableRange lvalue{{{0, 1}, {1, 2}}};
    EXPECT_EQ((CSRGraph<node, Unweighted>::fromEdges(3, lvalue).numberOfEdges()), 2);
    EXPECT_EQ((CSRGraph<node, Unweighted>::fromEdges(3, MutableRange{{{0, 1}}}).numberOfEdges()),
              1);
}

TEST(CSRGraphFromGraphGTest, testWeightedGraphToUnweighted) {
    // The weights of a weighted Graph are dropped when converting to an unweighted CSRGraph.
    Graph G(3, /*weighted=*/true);
    G.addEdge(0, 1, 2.5);
    G.addEdge(2, 1, 4.);
    const auto csr = CSRGraph<node, Unweighted>::fromGraph(G);
    static_assert(!decltype(csr)::isWeighted());
    EXPECT_EQ(csr.numberOfEdges(), 2);
    EXPECT_THAT(std::vector<node>(csr.neighborRange(1).begin(), csr.neighborRange(1).end()),
                ElementsAre(node{0}, node{2}));
}

TEST(CSRGraphFromGraphGTest, testUnweightedGraphToWeighted) {
    // An unweighted Graph reports weight 1 for every edge, and so does the converted graph.
    Graph G(3, /*weighted=*/false, /*directed=*/true);
    G.addEdge(0, 2);
    G.addEdge(0, 1);
    const auto csr = CSRGraph<node, edgeweight>::fromGraph(G);
    EXPECT_TRUE(csr.isDirected());
    EXPECT_THAT(weightedNeighborhoods(csr)[0], ElementsAre(Pair(node{1}, 1.), Pair(node{2}, 1.)));
}

TEST(CSRGraphFromGraphGTest, testUnweightedSourceTypeToWeighted) {
    // Same as above, for a source graph that is unweighted by type rather than at runtime.
    AdjListGraph<int, Unweighted> adj(3, /*weighted=*/false, /*directed=*/true);
    adj.addEdge(0, 2);
    adj.addEdge(0, 1);
    const auto csr = CSRGraph<int, float>::fromGraph(adj);
    EXPECT_THAT(weightedNeighborhoods(csr)[0], ElementsAre(Pair(1, 1.), Pair(2, 1.)));
}

TEST(CSRGraphFromGraphGTest, testFromCSRGraph) {
    // CSRGraphs can be converted to other node and weight types.
    const std::vector<WeightedEdgeT<node, double>> edgeList{{0, 1, 2.}, {2, 2, 3.}, {2, 0, 4.}};
    const auto source = CSRGraph<node, edgeweight>::fromEdges(3, edgeList);

    const auto weighted = CSRGraph<int, float>::fromGraph(source);
    EXPECT_EQ(weighted.numberOfEdges(), source.numberOfEdges());
    EXPECT_EQ(weighted.numberOfSelfLoops(), source.numberOfSelfLoops());
    EXPECT_THAT(weightedNeighborhoods(weighted)[2], ElementsAre(Pair(0, 4.), Pair(2, 3.)));

    const auto unweighted = CSRGraph<int, Unweighted>::fromGraph(source);
    EXPECT_EQ(unweighted.numberOfEdges(), source.numberOfEdges());
    EXPECT_EQ(unweighted.degree(0), source.degree(0));

    // An unweighted source gives every edge weight 1.
    EXPECT_THAT(weightedNeighborhoods(CSRGraph<node, edgeweight>::fromGraph(unweighted))[1],
                ElementsAre(Pair(node{0}, 1.)));
}

TEST(CSRGraphFromGraphGTest, testDeletedNodesBecomeIsolated) {
    // Node ids are kept, so a deleted node becomes an isolated node of the CSRGraph.
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(1, 3);
    G.addEdge(2, 3);
    G.removeEdge(2, 3);
    G.removeNode(2);
    ASSERT_EQ(G.numberOfNodes(), 3);

    const auto csr = CSRGraph<node, edgeweight>::fromGraph(G);
    EXPECT_EQ(csr.numberOfNodes(), G.upperNodeIdBound());
    EXPECT_EQ(csr.numberOfEdges(), G.numberOfEdges());
    EXPECT_EQ(csr.degree(2), 0);
    EXPECT_THAT(std::vector<node>(csr.neighborRange(1).begin(), csr.neighborRange(1).end()),
                ElementsAre(node{0}, node{3}));
}

} // namespace
} // namespace NetworKit

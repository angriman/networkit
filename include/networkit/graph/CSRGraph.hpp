#ifndef NETWORKIT_GRAPH_CSR_GRAPH_HPP_
#define NETWORKIT_GRAPH_CSR_GRAPH_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <omp.h>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/EdgeIterators.hpp>
#include <networkit/graph/EdgeLambda.hpp>
#include <networkit/graph/EdgeUtils.hpp>
#include <networkit/graph/GraphConcepts.hpp>
#include <networkit/graph/NodeIterators.hpp>

namespace NetworKit {

// Element of an edge list for a CSRGraph<NodeT, EdgeWeightT>: it has endpoints `u` and `v` (e.g.
// `EdgeT`), and a `weight` as well unless the graph is unweighted (e.g. `WeightedEdgeT`).
template <class EdgeType, class NodeT, class EdgeWeightT>
concept CSRGraphInputEdge = requires(const EdgeType &edge) {
    { edge.u } -> std::convertible_to<NodeT>;
    { edge.v } -> std::convertible_to<NodeT>;
} && (std::same_as<EdgeWeightT, Unweighted> || requires(const EdgeType &edge) {
                                { edge.weight } -> std::convertible_to<EdgeWeightT>;
                            });

// Graph in compressed sparse row format. The outgoing edges of a node are stored
// contiguously as interleaved `Edge{target, weight}` entries (so that traversals touch one cache
// line per neighbor), and each node's adjacency segment is sorted by target id. Graphs are meant to
// be built in bulk with fromEdges() or fromGraph().
template <GraphNode NodeType, GraphEdgeWeight EdgeWeightType>
class CSRGraph {
public:
    using NodeT = NodeType;
    using EdgeWeightT = EdgeWeightType;

    // For support of API: NetworKit::CSRGraph::NodeIterator
    using NodeIterator = NodeIteratorBase<CSRGraph, NodeT, EdgeWeightT>;
    using NodeRange = NodeRangeBase<CSRGraph, NodeT, EdgeWeightT>;

    // Iterators/ranges over the edges of the graph. `EdgeRange` yields plain edges (u, v), while
    // `EdgeWeightRange` also carries the edge weight. Both reuse the shared iterator machinery.
    using EdgeIterator = EdgeWeightTIterator<CSRGraph, NodeT, EdgeWeightT, EdgeT<NodeT>>;
    using EdgeWeightIterator =
        EdgeWeightTIterator<CSRGraph, NodeT, EdgeWeightT, WeightedEdgeT<NodeT, EdgeWeightT>>;
    using EdgeRange = EdgeWeightTRange<CSRGraph, NodeT, EdgeWeightT, EdgeT<NodeT>>;
    using EdgeWeightRange =
        EdgeWeightTRange<CSRGraph, NodeT, EdgeWeightT, WeightedEdgeT<NodeT, EdgeWeightT>>;

    // There is deliberately no constructor taking only `isDirected`: `CSRGraph(n)` would then be
    // ambiguous whenever `n` is an integer type other than `count` (e.g. a literal).
    explicit CSRGraph(count numNodes, bool isDirected = false) : directed(isDirected) {
        offset.resize(numNodes + 1);
    }

    CSRGraph() : CSRGraph(/*numNodes=*/0, /*isDirected=*/false) {}

    // Build a graph with `numNodes` nodes from a range of edges in O(n + m) (plus sorting each
    // adjacency segment); see CSRGraphInputEdge for the elements of the range. Weights are ignored
    // if the graph is unweighted. For undirected graphs each edge is given once, in either
    // direction. Parallel edges and self-loops are kept as given; parallel edges keep their
    // relative input order. The edge list is traversed twice, so it must be a forward range.
    // Defined here since old versions of clang cannot match an out-of-line definition to a
    // constraint that depends on the template parameters of the class.
    template <std::ranges::forward_range EdgeList>
        requires CSRGraphInputEdge<std::ranges::range_reference_t<EdgeList>, NodeT, EdgeWeightT>
    [[nodiscard]] static CSRGraph fromEdges(count numNodes, EdgeList &&edgeList,
                                            bool isDirected = false) {
        return build(numNodes, isDirected, [&edgeList](auto &&handle) -> void {
            for (const auto &edge : edgeList) {
                if constexpr (std::same_as<EdgeWeightT, Unweighted>) {
                    handle(static_cast<NodeT>(edge.u), static_cast<NodeT>(edge.v),
                           EdgeWeightInput{1});
                } else {
                    handle(static_cast<NodeT>(edge.u), static_cast<NodeT>(edge.v),
                           static_cast<EdgeWeightT>(edge.weight));
                }
            }
        });
    }

    // Build a graph with the same nodes and edges as `G`, an AdjListGraph, a CSRGraph or any graph
    // with the same iteration interface. Node ids are kept, so nodes of `G` that were deleted
    // become isolated nodes. The weights of `G` are dropped if this graph type is `Unweighted`, and
    // every edge gets weight 1 if `GraphT` is.
    template <class GraphT>
    [[nodiscard]] static CSRGraph fromGraph(const GraphT &G);

    count numberOfNodes() const noexcept {
        assert(!offset.empty());
        return offset.size() - 1;
    }
    count numberOfEdges() const noexcept {
        // Self-loops are stored once, every other undirected edge twice.
        assert(directed || (edges.size() - storedNumberOfSelfLoops) % 2 == 0);
        return directed ? edges.size() : (edges.size() + storedNumberOfSelfLoops) / 2;
    }
    count numberOfSelfLoops() const noexcept { return storedNumberOfSelfLoops; }
    NodeT upperNodeIdBound() const noexcept { return static_cast<NodeT>(numberOfNodes()); }
    bool isDirected() const noexcept { return directed; }
    static constexpr bool isWeighted() noexcept { return !std::same_as<EdgeWeightT, Unweighted>; }
    bool hasNode(NodeT u) const noexcept;

    // Range of the nodes of the graph, to be used in range-based for loops.
    NodeRange nodeRange() const noexcept { return NodeRange(*this); }

    struct Edge {
        NodeT target;
        // Optimize weights away in case of unweighted graphs (EdgeWeightT is Unweighted).
        [[no_unique_address]] EdgeWeightT weight;
    };

    // Number of (outgoing) neighbors of `u`, i.e., the size of its adjacency segment.
    count degree(NodeT u) const {
        assert(hasNode(u));
        return offset[u + 1] - offset[u];
    }

    // The `i`-th (outgoing) neighbor of `u`. `i` must be in [0, degree(u)); only checked by
    // assertions.
    NodeT getIthNeighbor(Unsafe, NodeT u, index i) const {
        assert(hasNode(u));
        assert(i < degree(u));
        return edges[offset[u] + i].target;
    }

    // Weight of the edge to the `i`-th (outgoing) neighbor of `u`. For unweighted graphs this
    // returns the (empty) `Unweighted` value. `i` must be in [0, degree(u)); only checked by
    // assertions.
    EdgeWeightT getIthNeighborWeight(Unsafe, NodeT u, index i) const {
        assert(hasNode(u));
        assert(i < degree(u));
        return edges[offset[u] + i].weight;
    }

    // Range over the edges of the graph. For undirected graphs each edge (u, v) is reported once,
    // with u <= v.
    EdgeRange edgeRange() const noexcept { return EdgeRange(*this); }

    // Range over the edges of the graph together with their weights (see edgeRange()).
    EdgeWeightRange edgeWeightRange() const noexcept { return EdgeWeightRange(*this); }

    // Iterator over the (outgoing) neighbors of a node, yielding the neighbor node ids. CSR stores
    // neighbors as `Edge`s, so this custom iterator projects each `Edge` onto its `target`.
    class NeighborIterator {
        const Edge *edgeIter = nullptr;

    public:
        using value_type = NodeT;
        // operator*() returns the projected target by value.
        using reference = value_type;
        using pointer = void;
        using iterator_category = std::forward_iterator_tag;
        using difference_type = ptrdiff_t;
        using self = NeighborIterator;

        NeighborIterator(const Edge *edgeIter) : edgeIter(edgeIter) {}
        // Required by the Python bindings.
        NeighborIterator() = default;

        NeighborIterator &operator++() {
            ++edgeIter;
            return *this;
        }

        NeighborIterator operator++(int) {
            const auto tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const NeighborIterator &rhs) const = default;

        NodeT operator*() const { return edgeIter->target; }
    };

    // Iterator over the (outgoing) neighbors of a node, yielding (neighbor, weight) pairs.
    class NeighborWeightIterator {
        const Edge *edgeIter = nullptr;

    public:
        using value_type = std::pair<NodeT, EdgeWeightT>;
        // operator*() builds the pair by value.
        using reference = value_type;
        using pointer = void;
        using iterator_category = std::forward_iterator_tag;
        using difference_type = ptrdiff_t;
        using self = NeighborWeightIterator;

        NeighborWeightIterator(const Edge *edgeIter) : edgeIter(edgeIter) {}
        // Required by the Python bindings.
        NeighborWeightIterator() = default;

        NeighborWeightIterator &operator++() {
            ++edgeIter;
            return *this;
        }

        NeighborWeightIterator operator++(int) {
            const auto tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const NeighborWeightIterator &rhs) const = default;

        std::pair<NodeT, EdgeWeightT> operator*() const {
            return {edgeIter->target, edgeIter->weight};
        }
    };

    // Range over the (outgoing) neighbors of `u`, to be used in range-based for loops.
    class NeighborRange {
        const CSRGraph *G = nullptr;
        NodeT u{};

    public:
        NeighborRange(const CSRGraph &G, NodeT u) : G(&G), u(u) { assert(G.hasNode(u)); }
        NeighborRange() {}

        NeighborIterator begin() const {
            assert(G != nullptr);
            return NeighborIterator(G->edges.data() + G->offset[u]);
        }
        NeighborIterator end() const {
            assert(G != nullptr);
            return NeighborIterator(G->edges.data() + G->offset[u + 1]);
        }
    };

    // Range over the (outgoing) neighbors of `u` together with the edge weights.
    class NeighborWeightRange {
        const CSRGraph *G = nullptr;
        NodeT u{};

    public:
        NeighborWeightRange(const CSRGraph &G, NodeT u) : G(&G), u(u) { assert(G.hasNode(u)); }
        NeighborWeightRange() {}

        NeighborWeightIterator begin() const {
            assert(G != nullptr);
            return NeighborWeightIterator(G->edges.data() + G->offset[u]);
        }
        NeighborWeightIterator end() const {
            assert(G != nullptr);
            return NeighborWeightIterator(G->edges.data() + G->offset[u + 1]);
        }
    };

    // Range over the (outgoing) neighbors of `u`.
    NeighborRange neighborRange(NodeT u) const { return NeighborRange(*this, u); }

    // Range over the (outgoing) neighbors of `u` together with the edge weights.
    NeighborWeightRange weightNeighborRange(NodeT u) const { return NeighborWeightRange(*this, u); }

    /* NODE ITERATORS */

    // Iterate over all nodes of the graph and call `handle(node)`.
    template <class Lambda>
    void forNodes(Lambda handle) const;

    // Iterate in parallel over all nodes of the graph and call `handle(node)`.
    template <class Lambda>
    void parallelForNodes(Lambda handle) const;

    /* EDGE ITERATORS */

    // The edge and neighborhood iterators below accept the same handle signatures as
    // AdjListGraph, except for those taking an edge id (see callEdgeLambda()):
    // `handle(u, v, weight)`, `handle(v, weight)`, `handle(u, v)` or `handle(v)`.

    // Iterate over all edges of the graph and call `handle` on each one. For undirected graphs
    // each edge is reported once, with u >= v, as in AdjListGraph::forEdges(). Note that
    // edgeRange() reports undirected edges with u <= v instead, again as in AdjListGraph.
    template <class Lambda>
    void forEdges(Lambda handle) const;

    // Iterate in parallel over all edges of the graph; see forEdges().
    template <class Lambda>
    void parallelForEdges(Lambda handle) const;

    /* NEIGHBORHOOD ITERATORS */

    // Iterate over the (outgoing) neighbors `v` of `u` and call `handle` on each edge (u, v). For
    // directed graphs only outgoing edges of `u` are considered. Same as forEdgesOf().
    template <class Lambda>
    void forNeighborsOf(NodeT u, Lambda handle) const;

    // Iterate over the incident (outgoing) edges (u, v) of `u` and call `handle` on each one.
    template <class Lambda>
    void forEdgesOf(NodeT u, Lambda handle) const;

    // Reserve memory for `numEdges` entries in the edge array (i.e. `numEdges` directed edges).
    void reserveEdges(size_t numEdges) { edges.reserve(numEdges); }

    // Type used to pass edge weights to `addEdge`: the graph's own weight type `EdgeWeightT` for
    // weighted graphs, and `edgeweight` for unweighted graphs so that a weight can still be
    // supplied (and then ignored).
    using EdgeWeightInput =
        std::conditional_t<std::same_as<EdgeWeightT, Unweighted>, edgeweight, EdgeWeightT>;

    // Slow helper (O(n + m) per call): prefer fromEdges() or fromGraph() for building graphs.
    // Insert an edge between nodes `u` and `v`. An edge weight can optionally be provided; for
    // unweighted graphs the weight is accepted but ignored (mirroring AdjListGraph). Parallel edges
    // are kept, after the existing edges with the same endpoints. Returns whether the edge was
    // added; currently always `true`, but a `bool` is returned to mirror AdjListGraph and to leave
    // room for a future multi-edge check.
    bool addEdge(NodeT u, NodeT v, EdgeWeightInput weight = EdgeWeightInput{1});

private:
    // Insert `v` (with weight `weight`) into the adjacency segment of `u` at its sorted position
    // (after existing entries with the same target), shifting the segments of all following nodes
    // one slot to the right so that the CSR layout stays contiguous. `weight` is ignored for
    // unweighted graphs.
    void addEdgeOneDirection(NodeT u, NodeT v, EdgeWeightInput weight);

    // Build an `Edge`; `weight` is ignored for unweighted graphs.
    static Edge makeEdge(NodeT target, [[maybe_unused]] EdgeWeightInput weight) {
        Edge edge;
        edge.target = target;
        if constexpr (!std::same_as<EdgeWeightT, Unweighted>) {
            edge.weight = weight;
        }
        return edge;
    }

    // Counting-sort construction shared by fromEdges() and fromGraph(). `forEachEdge(handle)` must
    // call `handle(u, v, weight)` once per input edge, and is invoked twice: first to count the
    // degrees, then to fill the adjacency segments.
    template <class ForEachEdge>
    static CSRGraph build(count numNodes, bool isDirected, ForEachEdge forEachEdge);

    // Invoke an edge handle with the signature it asks for. CSRGraph has no edge ids, so handles
    // taking one are rejected.
    template <class L>
    static decltype(auto) callEdgeHandle(L &handle, NodeT u, NodeT v, EdgeWeightT weight) {
        return callEdgeLambda<NodeT, EdgeWeightT, /*WithEdgeIds=*/false>(handle, u, v, weight,
                                                                         nullEdgeId);
    }

    // Report `u`'s contribution to a whole-graph edge iteration: each incident edge, but for
    // undirected graphs only once (with u >= v). Shared by forEdges() and parallelForEdges().
    template <class L>
    void forEdgesOfNode(NodeT u, L &handle) const;

protected:
    bool directed;
    count storedNumberOfSelfLoops = 0;
    std::vector<size_t> offset;
    std::vector<Edge> edges;
};

} // namespace NetworKit

#include <networkit/graph/CSRGraphImpl.hpp>

#endif // NETWORKIT_GRAPH_CSR_GRAPH_HPP_

#ifndef NETWORKIT_GRAPH_CSR_GRAPH_HPP_
#define NETWORKIT_GRAPH_CSR_GRAPH_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <limits>
#include <omp.h>
#include <type_traits>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/EdgeIterators.hpp>
#include <networkit/graph/EdgeUtils.hpp>
#include <networkit/graph/GraphConcepts.hpp>
#include <networkit/graph/NodeIterators.hpp>

namespace NetworKit {

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
class CSRGraph {
public:
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

    explicit CSRGraph(size_t numNodes, bool isDirected = false) : directed(isDirected) {
        offset.resize(numNodes + 1);
    }

    explicit CSRGraph(bool isDirected) : CSRGraph(/*numNodes=*/0, isDirected) {}

    CSRGraph() : CSRGraph(/*numNodes=*/0, /*isDirected=*/false) {}

    [[nodiscard]] count numberOfNodes() const noexcept {
        assert(!offset.empty());
        return offset.size() - 1;
    }
    [[nodiscard]] count numberOfEdges() const noexcept {
        // Self-loops are stored once, every other undirected edge twice.
        assert(directed || (edges.size() - storedNumberOfSelfLoops) % 2 == 0);
        return directed ? edges.size() : (edges.size() + storedNumberOfSelfLoops) / 2;
    }
    [[nodiscard]] count numberOfSelfLoops() const noexcept { return storedNumberOfSelfLoops; }
    [[nodiscard]] NodeT upperNodeIdBound() const noexcept {
        return static_cast<NodeT>(numberOfNodes());
    }
    [[nodiscard]] bool isDirected() const noexcept { return directed; }
    [[nodiscard]] static constexpr bool isWeighted() noexcept {
        return !std::same_as<EdgeWeightT, Unweighted>;
    }
    [[nodiscard]] bool hasNode(NodeT u) const noexcept;

    // Range of the nodes of the graph, to be used in range-based for loops.
    [[nodiscard]] NodeRange nodeRange() const noexcept { return NodeRange(*this); }

    struct Edge {
        NodeT target;
        // Optimize weights away in case of unweighted graphs (EdgeWeightT is Unweighted).
        [[no_unique_address]] EdgeWeightT weight;
    };

    // Number of (outgoing) neighbors of `u`, i.e. the size of its adjacency segment.
    [[nodiscard]] count degree(NodeT u) const {
        assert(hasNode(u));
        return offset[u + 1] - offset[u];
    }

    // The `i`-th (outgoing) neighbor of `u`. `i` must be in [0, degree(u)); not checked.
    [[nodiscard]] NodeT getIthNeighbor(Unsafe, NodeT u, index i) const {
        return edges[offset[u] + i].target;
    }

    // Weight of the edge to the `i`-th (outgoing) neighbor of `u`. For unweighted graphs this
    // returns the (empty) `Unweighted` value. `i` must be in [0, degree(u)); not checked.
    [[nodiscard]] EdgeWeightT getIthNeighborWeight(Unsafe, NodeT u, index i) const {
        return edges[offset[u] + i].weight;
    }

    // Range over the edges of the graph. For undirected graphs each edge (u, v) is reported once,
    // with u <= v.
    [[nodiscard]] EdgeRange edgeRange() const noexcept { return EdgeRange(*this); }

    // Range over the edges of the graph together with their weights (see edgeRange()).
    [[nodiscard]] EdgeWeightRange edgeWeightRange() const noexcept {
        return EdgeWeightRange(*this);
    }

    // Iterator over the (outgoing) neighbors of a node, yielding the neighbor node ids. CSR stores
    // neighbors as `Edge`s, so this custom iterator projects each `Edge` onto its `target`.
    class NeighborIterator {
        typename std::vector<Edge>::const_iterator edgeIter;

    public:
        using value_type = NodeT;
        using reference = value_type &;
        using pointer = value_type *;
        using iterator_category = std::forward_iterator_tag;
        using difference_type = ptrdiff_t;
        using self = NeighborIterator;

        NeighborIterator(typename std::vector<Edge>::const_iterator edgeIter)
            : edgeIter(edgeIter) {}
        // Required by the Python bindings; leaves the iterator uninitialized.
        NeighborIterator() {}

        NeighborIterator &operator++() {
            ++edgeIter;
            return *this;
        }

        NeighborIterator operator++(int) {
            const auto tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const NeighborIterator &rhs) const { return edgeIter == rhs.edgeIter; }

        bool operator!=(const NeighborIterator &rhs) const { return !(*this == rhs); }

        NodeT operator*() const { return edgeIter->target; }
    };

    // Iterator over the (outgoing) neighbors of a node, yielding (neighbor, weight) pairs.
    class NeighborWeightIterator {
        typename std::vector<Edge>::const_iterator edgeIter;

    public:
        using value_type = std::pair<NodeT, EdgeWeightT>;
        using reference = value_type &;
        using pointer = value_type *;
        using iterator_category = std::forward_iterator_tag;
        using difference_type = ptrdiff_t;
        using self = NeighborWeightIterator;

        NeighborWeightIterator(typename std::vector<Edge>::const_iterator edgeIter)
            : edgeIter(edgeIter) {}
        // Required by the Python bindings; leaves the iterator uninitialized.
        NeighborWeightIterator() {}

        NeighborWeightIterator &operator++() {
            ++edgeIter;
            return *this;
        }

        NeighborWeightIterator operator++(int) {
            const auto tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const NeighborWeightIterator &rhs) const {
            return edgeIter == rhs.edgeIter;
        }

        bool operator!=(const NeighborWeightIterator &rhs) const { return !(*this == rhs); }

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
            return NeighborIterator(G->edges.begin() + G->offset[u]);
        }
        NeighborIterator end() const {
            assert(G != nullptr);
            return NeighborIterator(G->edges.begin() + G->offset[u + 1]);
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
            return NeighborWeightIterator(G->edges.begin() + G->offset[u]);
        }
        NeighborWeightIterator end() const {
            assert(G != nullptr);
            return NeighborWeightIterator(G->edges.begin() + G->offset[u + 1]);
        }
    };

    // Range over the (outgoing) neighbors of `u`.
    [[nodiscard]] NeighborRange neighborRange(NodeT u) const { return NeighborRange(*this, u); }

    // Range over the (outgoing) neighbors of `u` together with the edge weights.
    [[nodiscard]] NeighborWeightRange weightNeighborRange(NodeT u) const {
        return NeighborWeightRange(*this, u);
    }

    /* NODE ITERATORS */

    // Iterate over all nodes of the graph and call `handle(node)`.
    template <typename Lambda>
    void forNodes(Lambda handle) const;

    // Iterate in parallel over all nodes of the graph and call `handle(node)`.
    template <typename Lambda>
    void parallelForNodes(Lambda handle) const;

    /* EDGE ITERATORS */

    // Iterate over all edges of the graph and call `handle(u, v)` or `handle(u, v, weight)`. For
    // undirected graphs each edge is reported once, with u <= v (see edgeRange()).
    template <typename Lambda>
    void forEdges(Lambda handle) const;

    // Iterate in parallel over all edges of the graph; see forEdges().
    template <typename Lambda>
    void parallelForEdges(Lambda handle) const;

    /* NEIGHBORHOOD ITERATORS */

    // Iterate over the (outgoing) neighbors of `u` and call `handle(v)` or `handle(v, weight)`. For
    // directed graphs only outgoing edges of `u` are considered.
    template <typename Lambda>
    void forNeighborsOf(NodeT u, Lambda handle) const;

    // Iterate over the incident (outgoing) edges of `u` and call `handle(u, v)` or
    // `handle(u, v, weight)`, where `v` is a neighbor of `u`.
    template <typename Lambda>
    void forEdgesOf(NodeT u, Lambda handle) const;

    // Reserve memory for `numEdges` entries in the edge array (i.e. `numEdges` directed edges).
    void reserveEdges(size_t numEdges) { edges.reserve(numEdges); }

    // Type used to pass edge weights to `addEdge`: the graph's own weight type `EdgeWeightT` for
    // weighted graphs, and `edgeweight` for unweighted graphs so that a weight can still be
    // supplied (and then ignored).
    using EdgeWeightInput =
        std::conditional_t<std::same_as<EdgeWeightT, Unweighted>, edgeweight, EdgeWeightT>;

    // Insert an edge between nodes `u` and `v`. An edge weight can optionally be provided; for
    // unweighted graphs the weight is accepted but ignored (mirroring AdjListGraph). Multi-edges
    // are not supported: inserting an already existing edge leaves the graph in an inconsistent
    // state. Returns whether the edge was added; currently always `true`, but a `bool` is returned
    // to mirror AdjListGraph and to leave room for a future multi-edge check.
    bool addEdge(NodeT u, NodeT v, EdgeWeightInput weight = EdgeWeightInput{1});

private:
    // Append `v` (with weight `weight`) to the adjacency segment of `u`, shifting the segments of
    // all following nodes one slot to the right so that the CSR layout stays contiguous. `weight`
    // is ignored for unweighted graphs.
    void addEdgeOneDirection(NodeT u, NodeT v, EdgeWeightInput weight);

    // Invoke an edge handle, choosing the `(u, v)` or `(u, v, weight)` signature based on arity.
    template <typename L>
    void callEdgeHandle(L &handle, NodeT u, NodeT v, EdgeWeightT weight) const {
        if constexpr (std::is_invocable_v<L, NodeT, NodeT, EdgeWeightT>) {
            handle(u, v, weight);
        } else {
            static_assert(std::is_invocable_v<L, NodeT, NodeT>,
                          "Edge handle must be callable as (node, node) or (node, node, weight).");
            handle(u, v);
        }
    }

    // Invoke a neighbor handle, choosing the `(v)` or `(v, weight)` signature based on arity.
    template <typename L>
    void callNeighborHandle(L &handle, NodeT v, EdgeWeightT weight) const {
        if constexpr (std::is_invocable_v<L, NodeT, EdgeWeightT>) {
            handle(v, weight);
        } else {
            static_assert(std::is_invocable_v<L, NodeT>,
                          "Neighbor handle must be callable as (node) or (node, weight).");
            handle(v);
        }
    }

    // Report `u`'s contribution to a whole-graph edge iteration: each incident edge, but for
    // undirected graphs only once (with u <= v). Shared by forEdges() and parallelForEdges().
    template <typename L>
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

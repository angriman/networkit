#ifndef NETWORKIT_GRAPH_CSR_GRAPH_IMPL_HPP_
#define NETWORKIT_GRAPH_CSR_GRAPH_IMPL_HPP_

namespace NetworKit {

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
bool CSRGraph<NodeT, EdgeWeightT>::hasNode(NodeT u) const noexcept {
    if constexpr (std::is_signed<NodeT>::value) {
        return u >= 0 && static_cast<count>(u) < numberOfNodes();
    } else {
        return static_cast<count>(u) < numberOfNodes();
    }
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
bool CSRGraph<NodeT, EdgeWeightT>::addEdge(NodeT u, NodeT v, EdgeWeightInput weight) {
    assert(hasNode(u) && hasNode(v));
    addEdgeOneDirection(u, v, weight);
    if (u == v) {
        // Self-loops are stored only once, even in undirected graphs.
        ++storedNumberOfSelfLoops;
    } else if (!directed) {
        addEdgeOneDirection(v, u, weight);
    }
    return true;
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
void CSRGraph<NodeT, EdgeWeightT>::addEdgeOneDirection(NodeT u, NodeT v,
                                                       [[maybe_unused]] EdgeWeightInput weight) {
    Edge newEdge;
    newEdge.target = v;
    // For weighted graphs `EdgeWeightInput` is `EdgeWeightT`, so the weight is stored directly.
    // For unweighted graphs (EdgeWeightT is Unweighted) the weight is accepted but not stored.
    if constexpr (!std::same_as<EdgeWeightT, Unweighted>) {
        newEdge.weight = weight;
    }

    // Insert the new edge at the end of u's adjacency segment. std::vector::insert shifts every
    // following entry one slot to the right, keeping the CSR layout contiguous.
    edges.insert(edges.begin() + offset[u + 1], newEdge);

    // Every node after u now starts one slot further along the edge array.
    std::for_each(offset.begin() + u + 1, offset.end(), [](size_t &o) { ++o; });
}

/* NODE ITERATORS */

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <typename Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forNodes(Lambda handle) const {
    const NodeT n = static_cast<NodeT>(numberOfNodes());
    for (NodeT u = 0; u < n; ++u) {
        handle(u);
    }
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <typename Lambda>
void CSRGraph<NodeT, EdgeWeightT>::parallelForNodes(Lambda handle) const {
    const auto n = static_cast<omp_index>(numberOfNodes());
#pragma omp parallel for
    for (omp_index u = 0; u < n; ++u) {
        handle(static_cast<NodeT>(u));
    }
}

/* EDGE ITERATORS */

// The neighborhood loops below hoist a pointer to `u`'s adjacency segment and its degree into
// locals before iterating. The handle is an opaque call, so the compiler cannot prove it leaves
// `edges`/`offset` unchanged; without hoisting it would reload those vector base pointers and
// recompute the degree on every iteration. Hoisting yields the same tight loop as the range API.
// This is also why we read `segment[i]` directly instead of calling getIthNeighbor(u, i) /
// getIthNeighborWeight(u, i): each of those recomputes `edges[offset[u] + i]` from scratch, so
// using them here would reintroduce the per-iteration `offset[u]` reloads we are avoiding.

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <typename Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forEdgesOfNode(NodeT u, Lambda &handle) const {
    const Edge *segment = edges.data() + offset[u];
    const count deg = degree(u);
    for (index i = 0; i < deg; ++i) {
        const NodeT v = segment[i].target;
        // Undirected edges are stored in both directions; report each once (u <= v).
        if (directed || u <= v) {
            callEdgeHandle(handle, u, v, segment[i].weight);
        }
    }
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <typename Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forEdges(Lambda handle) const {
    const NodeT n = static_cast<NodeT>(numberOfNodes());
    for (NodeT u = 0; u < n; ++u) {
        forEdgesOfNode(u, handle);
    }
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <typename Lambda>
void CSRGraph<NodeT, EdgeWeightT>::parallelForEdges(Lambda handle) const {
    const auto n = static_cast<omp_index>(numberOfNodes());
#pragma omp parallel for schedule(guided)
    for (omp_index u_ = 0; u_ < n; ++u_) {
        forEdgesOfNode(static_cast<NodeT>(u_), handle);
    }
}

/* NEIGHBORHOOD ITERATORS */

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <typename Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forNeighborsOf(NodeT u, Lambda handle) const {
    forEdgesOf(
        u, [&handle, this](NodeT, NodeT v, EdgeWeightT w) { callNeighborHandle(handle, v, w); });
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <typename Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forEdgesOf(NodeT u, Lambda handle) const {
    assert(hasNode(u));
    const Edge *segment = edges.data() + offset[u];
    const count deg = degree(u);
    for (index i = 0; i < deg; ++i) {
        callEdgeHandle(handle, u, segment[i].target, segment[i].weight);
    }
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_CSR_GRAPH_IMPL_HPP_

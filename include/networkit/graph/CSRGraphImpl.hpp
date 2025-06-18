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
void CSRGraph<NodeT, EdgeWeightT>::addEdgeOneDirection(NodeT u, NodeT v, EdgeWeightInput weight) {
    // Keep u's segment sorted by target: insert after the entries that are not greater than v.
    const auto pos =
        std::ranges::upper_bound(targets.begin() + offset[u], targets.begin() + offset[u + 1], v);
    // std::vector::insert shifts every following entry one slot to the right, keeping the CSR
    // layout contiguous.
    if constexpr (isWeighted()) {
        weights.insert(weights.begin() + (pos - targets.begin()), weight);
    }
    targets.insert(pos, v);

    // Every node after u now starts one slot further along the edge arrays.
    std::ranges::for_each(offset.begin() + u + 1, offset.end(), [](size_t &o) -> void { ++o; });
}

/* BULK CONSTRUCTION */

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <class ForEachEdge>
CSRGraph<NodeT, EdgeWeightT> CSRGraph<NodeT, EdgeWeightT>::build(count numNodes, bool isDirected,
                                                                 ForEachEdge forEachEdge) {
    CSRGraph G(numNodes, isDirected);

    // Count the entries of every adjacency segment, shifted by one so that the prefix sum below
    // turns the counts into offsets. Undirected edges are stored in both directions, except
    // self-loops, which are stored once.
    forEachEdge([&G, isDirected](NodeT u, NodeT v, EdgeWeightInput) -> void {
        assert(G.hasNode(u) && G.hasNode(v));
        ++G.offset[u + 1];
        if (u == v) {
            ++G.storedNumberOfSelfLoops;
        } else if (!isDirected) {
            ++G.offset[v + 1];
        }
    });
    std::partial_sum(G.offset.begin(), G.offset.end(), G.offset.begin());

    // Scatter the edges into their segments. `next[u]` is the next free slot of u's segment.
    G.targets.resize(G.offset.back());
    if constexpr (isWeighted()) {
        G.weights.resize(G.offset.back());
    }
    std::vector<size_t> next(G.offset.begin(), G.offset.end() - 1);
    forEachEdge([&G, &next, isDirected](NodeT u, NodeT v, EdgeWeightInput weight) -> void {
        G.setEntry(next[u]++, v, weight);
        if (u != v && !isDirected) {
            G.setEntry(next[v]++, u, weight);
        }
    });

    G.sortAdjacencySegments();
    return G;
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
void CSRGraph<NodeT, EdgeWeightT>::sortAdjacencySegments() {
    const auto n = static_cast<omp_index>(numberOfNodes());
    if constexpr (isWeighted()) {
        // Targets and weights must move together, so sort copies of the pairs of each segment.
#pragma omp parallel
        {
            std::vector<std::pair<NodeT, EdgeWeightT>> segment;
#pragma omp for schedule(guided)
            for (omp_index u = 0; u < n; ++u) {
                const size_t first = offset[u];
                const size_t last = offset[u + 1];
                segment.clear();
                for (size_t i = first; i < last; ++i) {
                    segment.emplace_back(targets[i], weights[i]);
                }
                std::ranges::stable_sort(segment, {}, &std::pair<NodeT, EdgeWeightT>::first);
                for (size_t i = first; i < last; ++i) {
                    std::tie(targets[i], weights[i]) = segment[i - first];
                }
            }
        }
    } else {
        // Parallel edges are indistinguishable here, so stability does not matter.
#pragma omp parallel for schedule(guided)
        for (omp_index u = 0; u < n; ++u) {
            std::ranges::sort(targets.begin() + offset[u], targets.begin() + offset[u + 1]);
        }
    }
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <class GraphT>
CSRGraph<NodeT, EdgeWeightT> CSRGraph<NodeT, EdgeWeightT>::fromGraph(const GraphT &G) {
    using SrcNodeT = typename GraphT::NodeT;
    using SrcWeightT = typename GraphT::EdgeWeightT;
    return build(G.upperNodeIdBound(), G.isDirected(), [&G](auto &&handle) -> void {
        if constexpr (std::same_as<EdgeWeightT, Unweighted>
                      || std::same_as<SrcWeightT, Unweighted>) {
            G.forEdges([&handle](SrcNodeT u, SrcNodeT v) -> void {
                handle(static_cast<NodeT>(u), static_cast<NodeT>(v), EdgeWeightInput{1});
            });
        } else {
            G.forEdges([&handle](SrcNodeT u, SrcNodeT v, SrcWeightT w) -> void {
                handle(static_cast<NodeT>(u), static_cast<NodeT>(v), static_cast<EdgeWeightT>(w));
            });
        }
    });
}

/* NODE ITERATORS */

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <class Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forNodes(Lambda handle) const {
    const NodeT n = static_cast<NodeT>(numberOfNodes());
    for (NodeT u = 0; u < n; ++u) {
        handle(u);
    }
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <class Lambda>
void CSRGraph<NodeT, EdgeWeightT>::parallelForNodes(Lambda handle) const {
    const auto n = static_cast<omp_index>(numberOfNodes());
#pragma omp parallel for
    for (omp_index u = 0; u < n; ++u) {
        handle(static_cast<NodeT>(u));
    }
}

/* EDGE ITERATORS */

// The neighborhood loops below hoist pointers to `u`'s adjacency segment and its degree into locals
// before iterating. The handle is an opaque call, so the compiler cannot prove it leaves
// `targets`/`weights`/`offset` unchanged; without hoisting it would reload those vector base
// pointers and recompute the degree on every iteration. Hoisting yields the same tight loop as the
// range API. This is also why we read `segment[i]` directly instead of calling getIthNeighbor(u, i)
// / getIthNeighborWeight(u, i): each of those recomputes the segment from `offset[u]`, so using
// them here would reintroduce the per-iteration reloads we are avoiding.

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <class Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forEdgesOfNode(NodeT u, Lambda &handle) const {
    const NodeT *segment = targets.data() + offset[u];
    const EdgeWeightT *segmentWeights = weightSegment(u);
    const count deg = degree(u);
    for (index i = 0; i < deg; ++i) {
        const NodeT v = segment[i];
        // Undirected edges are stored in both directions; report each once (u >= v), like
        // AdjListGraph does.
        if (directed || u >= v) {
            callEdgeHandle(handle, u, v, weightAt(segmentWeights, i));
        }
    }
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <class Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forEdges(Lambda handle) const {
    const NodeT n = static_cast<NodeT>(numberOfNodes());
    for (NodeT u = 0; u < n; ++u) {
        forEdgesOfNode(u, handle);
    }
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <class Lambda>
void CSRGraph<NodeT, EdgeWeightT>::parallelForEdges(Lambda handle) const {
    const auto n = static_cast<omp_index>(numberOfNodes());
#pragma omp parallel for schedule(guided)
    for (omp_index u_ = 0; u_ < n; ++u_) {
        forEdgesOfNode(static_cast<NodeT>(u_), handle);
    }
}

/* NEIGHBORHOOD ITERATORS */

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <class Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forNeighborsOf(NodeT u, Lambda handle) const {
    forEdgesOf(u, handle);
}

template <GraphNode NodeT, GraphEdgeWeight EdgeWeightT>
template <class Lambda>
void CSRGraph<NodeT, EdgeWeightT>::forEdgesOf(NodeT u, Lambda handle) const {
    assert(hasNode(u));
    const NodeT *segment = targets.data() + offset[u];
    const EdgeWeightT *segmentWeights = weightSegment(u);
    const count deg = degree(u);
    for (index i = 0; i < deg; ++i) {
        callEdgeHandle(handle, u, segment[i], weightAt(segmentWeights, i));
    }
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_CSR_GRAPH_IMPL_HPP_

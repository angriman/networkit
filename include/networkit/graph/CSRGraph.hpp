#ifndef NETWORKIT_GRAPH_CSR_GRAPH_HPP_
#define NETWORKIT_GRAPH_CSR_GRAPH_HPP_

#include <algorithm>
#include <cassert>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/GeneralGraph.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

template <class NodeType, class EdgeWeightType>
class CSRGraphBase : public GeneralGraph<NodeType, EdgeWeightType> {
    static_assert(std::is_integral<NodeType>::value, "Integral required for NodeType.");

protected:
    using NodeRange = NodeRangeBase<CSRGraphBase<NodeType, EdgeWeightType>, NodeType>;

    std::vector<size_t> offset;
    bool directed;

    // Returns true if `n` nodes can fit into an integer of type `NodeType`. Returns false
    // otherwise.
    static bool nodesFitInGraph(size_t n) {
        return n == 0 || static_cast<size_t>(std::numeric_limits<NodeType>::max()) >= n - 1;
    }

public:
    CSRGraphBase() : directed(false) { offset.resize(1); }

    CSRGraphBase(size_t numNodes, bool isDirected = false) : directed(isDirected) {
        assert(nodesFitInGraph(numNodes));
        offset.resize(numNodes + 1);
    }

    count numberOfNodes() const noexcept override { return this->offset.size() - 1; }
    index upperNodeIdBound() const noexcept override { return numberOfNodes(); }

    bool isDirected() const noexcept { return this->directed; }

    bool hasNode(NodeType u) const override {
        if constexpr (std::is_signed<NodeType>::value) {
            return u >= 0 && static_cast<count>(u) < numberOfNodes();
        } else {
            return static_cast<count>(u) < numberOfNodes();
        }
    }

    NodeType addNode() {
        assert(this->numberOfNodes() < std::numeric_limits<NodeType>::max());
        this->offset.push_back(this->offset.back());
        return this->offset.back();
    }

    template <class Lambda>
    void parallelForNodes(Lambda lambda) const;

    NodeRange nodeRange() const noexcept { return NodeRange(*this); }
};

// Forward declaration of the `CSRGraph` class to enable specialization below.
template <class NodeType, class EdgeWeightType = void>
class CSRGraph;

/*
 * `CSRGraph` class for weighted graphs.
 */
template <class NodeType, class EdgeWeightType>
class CSRGraph : public CSRGraphBase<NodeType, EdgeWeightType> {
    using NodeIterator = NodeIteratorBase<CSRGraph<NodeType, EdgeWeightType>, NodeType>;

    template <bool InEdges>
    class NeighborRange;
    using OutNeighborRange = NeighborRange<false>;
    class NeighborWeightIterator;
    class NeighborWeightRange;

    struct Edge {
        NodeType neighbor;
        EdgeWeightType weight;
    };

    std::vector<Edge> edges;

    friend class NeighborIteratorBase<NodeType, Edge>;

public:
    // Iterators to iterate over the neighbors of a given node.
    class NeighborIterator;

    // Use the same default constructor as `CSRGraphBase`.
    using CSRGraphBase<NodeType, EdgeWeightType>::CSRGraphBase;

    CSRGraph(size_t numNodes, bool isDirected = false)
        : CSRGraphBase<NodeType, EdgeWeightType>::CSRGraphBase(numNodes, isDirected) {}

    CSRGraph(const Graph &dynGraph);

    void reserveEdges(size_t numEdges) { edges.reserve(numEdges + 1); }

    count numberOfEdges() const noexcept override {
        return this->directed ? edges.size() : edges.size() / 2;
    }
    index upperEdgeIdBound() const noexcept override { return numberOfEdges(); }

    OutNeighborRange neighborRange(NodeType u) const { return OutNeighborRange(*this, u); }
    NeighborWeightRange neighborWeightRange(NodeType u) const {
        return NeighborWeightRange(*this, u);
    }

    void addEdge(NodeType u, NodeType v,
                 EdgeWeightType weight = EdgeWeightType{defaultEdgeWeight}) {
        assert(this->hasNode(u) && this->hasNode(v));
        addEdgeOneDirection(u, v, weight);
        if (!this->isDirected()) {
            addEdgeOneDirection(v, u, weight);
        }
    }

protected:
    void addEdgeOneDirection(NodeType u, NodeType v, EdgeWeightType weight);
};

/*
 * `CSRGraph` class specialized for unweighted graphs.
 */
template <class NodeType>
class CSRGraph<NodeType, void> : public CSRGraphBase<NodeType, void> {
    std::vector<NodeType> edges;

    friend class NeighborIteratorBase<NodeType, NodeType>;

protected:
    using NodeIterator = NodeIteratorBase<CSRGraph<NodeType, void>, NodeType>;
    using NodeRange = typename CSRGraphBase<NodeType, void>::NodeRange;

    template <bool InEdges>
    class NeighborRange;
    using OutNeighborRange = NeighborRange<false>;
    class NeighborWeightIterator;
    class NeighborWeightRange;

public:
    // Iterators to iterate over the neighbors of a given node.
    class NeighborIterator;

    // Use the same default constructor as `CSRGraphBase`.
    using CSRGraphBase<NodeType, void>::CSRGraphBase;

    CSRGraph(size_t numNodes, bool isDirected = false)
        : CSRGraphBase<NodeType, void>::CSRGraphBase(numNodes, isDirected) {}

    CSRGraph(const Graph &dynGraph);

    void reserveEdges(size_t numEdges) { edges.reserve(numEdges + 1); }

    count numberOfEdges() const noexcept override {
        return this->directed ? edges.size() : edges.size() / 2;
    }
    index upperEdgeIdBound() const noexcept override { return numberOfEdges(); }

    OutNeighborRange neighborRange(NodeType u) const { return OutNeighborRange(*this, u); }
    NeighborWeightRange neighborWeightRange(NodeType u) const {
        return NeighborWeightRange(*this, u);
    }

    void addEdge(NodeType u, NodeType v) {
        assert(this->hasNode(u) && this->hasNode(v));
        addEdgeOneDirection(u, v);
        if (!this->isDirected()) {
            addEdgeOneDirection(v, u);
        }
    }

    template <class EdgeWeightType>
    void addEdge(NodeType u, NodeType v, EdgeWeightType) {
        WARN("Adding weighted edge to unweighted graph, the edge weight will be ignored!");
        addEdge(u, v);
    }

protected:
    void addEdgeOneDirection(NodeType u, NodeType v);
};

template <class NodeType, class EdgeWeightType>
template <class Lambda>
void CSRGraphBase<NodeType, EdgeWeightType>::parallelForNodes(Lambda lambda) const {
    NodeType numNodes = this->numberOfNodes();
#pragma omp parallel for schedule(dynamic)
    for (NodeType u = 0; u < numNodes; ++u) {
        lambda(u);
    }
}

/**
 * Wrapper class to iterate over a range of the nodes of a graph.
 */

template <class NodeType>
class CSRGraph<NodeType, void>::NeighborIterator : public NeighborIteratorBase<NodeType, NodeType> {
public:
    // Own type.
    using self = NeighborIterator;

    // Using the same constructor as `NeighborIteratorBase`.
    using NeighborIteratorBase<NodeType, NodeType>::NeighborIteratorBase;

    NodeType operator*() const noexcept {
        assert(this->neighIter != this->edges->end());
        return *this->neighIter;
    }
};

template <class NodeType, class EdgeWeightType>
class CSRGraph<NodeType, EdgeWeightType>::NeighborIterator
    : public NeighborIteratorBase<NodeType, Edge> {
public:
    // Own type.
    using self = NeighborIterator;

    // Using the same constructor as `NeighborIteratorBase`.
    using NeighborIteratorBase<NodeType, Edge>::NeighborIteratorBase;

    NodeType operator*() const noexcept {
        assert(this->neighIter != this->edges->end());
        return (*this->neighIter).neighbor;
    }
};

template <class NodeType>
class CSRGraph<NodeType, void>::NeighborWeightIterator
    : public NeighborIteratorBase<NodeType, NodeType> {
public:
    // Own type
    using self = NeighborWeightIterator;

    // Using the same constructor as `NeighborIteratorBase`.
    using NeighborIteratorBase<NodeType, NodeType>::NeighborIteratorBase;

    std::pair<NodeType, edgeweight> operator*() const noexcept {
        assert(this->neighIter != this->edges->end());
        return {*this->neighIter, defaultEdgeWeight};
    }
};

template <class NodeType, class EdgeWeightType>
class CSRGraph<NodeType, EdgeWeightType>::NeighborWeightIterator
    : public NeighborIteratorBase<NodeType, Edge> {
public:
    // Own type
    using self = NeighborWeightIterator;

    // Using the same constructor as `NeighborIteratorBase`.
    using NeighborIteratorBase<NodeType, Edge>::NeighborIteratorBase;

    std::pair<NodeType, EdgeWeightType> operator*() const noexcept {
        assert(this->neighIter != this->edges->end());
        return {(*this->neighIter).neighbor, (*this->neighIter).weight};
    }
};

template <class NodeType>
template <bool InEdges>
class CSRGraph<NodeType, void>::NeighborRange
    : NeighborRangeBase<CSRGraph<NodeType, void>, NodeType, InEdges> {
    friend class NeighborIterator;

public:
    using NeighborRangeBase<CSRGraph<NodeType, void>, NodeType, InEdges>::NeighborRangeBase;

    NeighborIterator begin() const {
        assert(this->G);
        return {this->G->edges, this->G->offset[this->source]};
    }

    NeighborIterator end() const {
        assert(this->G);
        return {this->G->edges, this->G->offset[this->source + 1]};
    }
};

template <class NodeType, class EdgeWeightType>
template <bool InEdges>
class CSRGraph<NodeType, EdgeWeightType>::NeighborRange
    : public NeighborRangeBase<CSRGraph<NodeType, EdgeWeightType>, NodeType, InEdges> {
    friend class NeighborIterator;

public:
    using NeighborRangeBase<CSRGraph<NodeType, EdgeWeightType>, NodeType,
                            InEdges>::NeighborRangeBase;

    NeighborIterator begin() const override {
        assert(this->G);
        return {this->G->edges, this->G->offset[this->source]};
    }

    NeighborIterator end() const override {
        assert(this->G);
        return {this->G->edges, this->G->offset[this->source + 1]};
    }
};

template <class NodeType>
class CSRGraph<NodeType, void>::NeighborWeightRange {
    friend class NeighborWeightIterator;

    const CSRGraph<NodeType, void> *G = nullptr;
    NodeType u;

public:
    NeighborWeightRange(const CSRGraph<NodeType, void> &G, NodeType u) : G(&G), u(u) {
        assert(G.hasNode(u));
    }

    NeighborWeightIterator begin() const {
        assert(G != nullptr);
        return {G->edges, G->offset[u]};
    }

    NeighborWeightIterator end() const {
        assert(G);
        return {G->edges, G->offset[u + 1]};
    }
};

template <class NodeType, class EdgeWeightType>
class CSRGraph<NodeType, EdgeWeightType>::NeighborWeightRange {
    friend class NeighborWeightIterator;

    const CSRGraph<NodeType, EdgeWeightType> *G = nullptr;
    NodeType u;

public:
    NeighborWeightRange(const CSRGraph<NodeType, EdgeWeightType> &G, NodeType u) : G(&G), u(u) {
        assert(G.hasNode(u));
    }

    NeighborWeightIterator begin() const {
        assert(G);
        return {G->edges, G->offset[u]};
    }

    NeighborWeightIterator end() const {
        assert(G);
        return {G->edges, G->offset[u + 1]};
    }
};

template <class NodeType>
CSRGraph<NodeType, void>::CSRGraph(const Graph &dynGraph) {
    assert(dynGraph.numberOfNodes() == dynGraph.upperNodeIdBound()
           && "Error, the input graph contains non-existing nodes.");

    this->directed = dynGraph.isDirected();
    this->offset.resize(dynGraph.upperNodeIdBound() + 1);
    edges.resize(dynGraph.isDirected() ? dynGraph.numberOfEdges() : dynGraph.numberOfEdges() * 2);

    size_t degreeSum = 0;

    for (size_t i = 0; i < dynGraph.upperNodeIdBound(); ++i) {
        this->offset[i] = degreeSum;
        for (node neighbor : dynGraph.neighborRange(i)) {
            edges[degreeSum] = static_cast<NodeType>(neighbor);
            ++degreeSum;
        }
    }

    this->offset.back() = degreeSum;
}

// Put helper functions within the anonymous namespace.
namespace {

template <class NodeType, class EdgeType>
void addEdgeInEdgesVector(std::vector<size_t> &offset, std::vector<EdgeType> &edges, NodeType u,
                          EdgeType newEdge) {
    // Position in the `edges` vector where the new edge should be stored.
    const size_t edgeOffset = offset[u + 1];

    // Allocate space for the new edge.
    edges.resize(edges.size() + 1);

    // Unless the edge has to be inserted in the back of `edges`, left rotate
    // all the edges by 1.
    if (edgeOffset < edges.size() - 1) {
        std::rotate(edges.rbegin(), edges.rbegin() + 1,
                    edges.rbegin() + edges.size() - 1 - edgeOffset);
    }

    // Store the new edge and update offsets of the next nodes.
    edges[edgeOffset] = newEdge;
    std::for_each(offset.begin() + u + 1, offset.end(), [](size_t &offsetValue) { ++offsetValue; });
}
} // namespace

template <class NodeType>
void CSRGraph<NodeType, void>::addEdgeOneDirection(NodeType u, NodeType v) {
    addEdgeInEdgesVector(this->offset, this->edges, u, v);
}

template <class NodeType, class EdgeWeightType>
void CSRGraph<NodeType, EdgeWeightType>::addEdgeOneDirection(NodeType u, NodeType v,
                                                             EdgeWeightType weight) {
    Edge newEdge;
    newEdge.neighbor = v;
    newEdge.weight = weight;
    addEdgeInEdgesVector(this->offset, this->edges, u, newEdge);
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_CSR_GRAPH_HPP_

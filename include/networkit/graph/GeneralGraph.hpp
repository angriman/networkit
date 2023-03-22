#ifndef NETWORKIT_GRAPH_GENERAL_GRAPH_HPP_
#define NETWORKIT_GRAPH_GENERAL_GRAPH_HPP_

#include <networkit/Globals.hpp>

namespace NetworKit {

template <class NodeType = node, class EdgeWeightType = edgeweight>
class GeneralGraph {
public:
    virtual count numberOfNodes() const = 0;
    virtual index upperNodeIdBound() const = 0;
    virtual count numberOfEdges() const = 0;
    virtual index upperEdgeIdBound() const = 0;
    virtual bool hasNode(NodeType u) const = 0;
};

template <class GraphType, class NodeType>
class NodeIteratorBase {
    const GraphType *G;
    NodeType u;

public:
    // The value type of the nodes (i.e. nodes). Returned by
    // operator*().
    using value_type = NodeType;

    // Reference to the value_type, required by STL.
    using reference = value_type &;

    // Pointer to the value_type, required by STL.
    using pointer = value_type *;

    // STL iterator category.
    using iterator_category = std::forward_iterator_tag;

    // Signed integer type of the result of subtracting two pointers,
    // required by STL.
    using difference_type = ptrdiff_t;

    // Own type.
    using self = NodeIteratorBase;

    NodeIteratorBase() : G(nullptr) {}

    NodeIteratorBase(const GraphType *G, NodeType u) : G(G), u(u) {
        if (!G->hasNode(u) && static_cast<index>(u) < G->upperNodeIdBound()) {
            ++(*this);
        }
    }

    ~NodeIteratorBase() = default;

    NodeIteratorBase &operator++() {
        assert(static_cast<index>(u) < G->upperNodeIdBound());
        do {
            ++u;
        } while (!(G->hasNode(u) || static_cast<index>(u) >= G->upperNodeIdBound()));
        return *this;
    }

    NodeIteratorBase operator++(int) {
        const auto tmp = *this;
        ++(*this);
        return tmp;
    }

    NodeIteratorBase operator--() {
        assert(u > 0);
        do {
            --u;
        } while (!G->hasNode(u));
        return *this;
    }

    NodeIteratorBase operator--(int) {
        const auto tmp = *this;
        --(*this);
        return tmp;
    }

    bool operator==(const NodeIteratorBase &rhs) const noexcept { return u == rhs.u; }

    bool operator!=(const NodeIteratorBase &rhs) const noexcept { return !(*this == rhs); }

    NodeType operator*() const noexcept {
        assert(static_cast<index>(u) < G->upperNodeIdBound());
        return u;
    }
};

/**
 * Wrapper class to iterate over a range of the nodes of a graph.
 */
template <class GraphType, class NodeType>
class NodeRangeBase {
    const GraphType *G;

public:
    NodeRangeBase() : G(nullptr) {}
    NodeRangeBase(const GraphType &G) : G(&G) {}
    ~NodeRangeBase() = default;

    NodeIteratorBase<GraphType, NodeType> begin() const noexcept {
        assert(G != nullptr);
        return NodeIteratorBase<GraphType, NodeType>(G, NodeType{0});
    }

    NodeIteratorBase<GraphType, NodeType> end() const noexcept {
        assert(G != nullptr);
        return NodeIteratorBase<GraphType, NodeType>(G, G->upperNodeIdBound());
    }
};

template <class NodeType>
class NeighborIterator;

template <class NodeType, class EdgeType>
class NeighborIteratorBase {
    friend class NeighborIterator<NodeType>;

protected:
    const std::vector<EdgeType> *edges;

public:
    typename std::vector<EdgeType>::const_iterator neighIter;

    // The value type of the nodes. Returned by operator*().
    using value_type = NodeType;

    // Reference to the value_type, required by STL.
    using reference = value_type &;

    // Pointer to the value_type, required by STL.
    using pointer = value_type *;

    // STL iterator category.
    using iterator_category = std::forward_iterator_tag;

    // Signed integer type of the result of subtracting two pointers,
    // required by STL.
    using difference_type = ptrdiff_t;

    NeighborIteratorBase(typename std::vector<EdgeType>::const_iterator neighIter)
        : edges(nullptr), neighIter(neighIter) {}

    NeighborIteratorBase(const std::vector<EdgeType> &edges, size_t offset = 0)
        : edges(&edges), neighIter(edges.begin() + offset) {}

    ~NeighborIteratorBase() = default;

    NeighborIteratorBase &operator++() {
        if (edges != nullptr) {
            assert(neighIter != edges->end());
        }
        ++neighIter;
        return *this;
    }

    NeighborIteratorBase operator++(int) {
        const auto tmp = *this;
        ++(*this);
        return tmp;
    }

    NeighborIteratorBase operator--() {
        if (edges != nullptr) {
            assert(neighIter != edges->begin());
        }
        --neighIter;
        return *this;
    }

    NeighborIteratorBase operator--(int) {
        const auto tmp = *this;
        --(*this);
        return tmp;
    }

    bool operator==(const NeighborIteratorBase &rhs) const noexcept {
        return neighIter == rhs.neighIter;
    }

    bool operator!=(const NeighborIteratorBase &rhs) const noexcept {
        return !(neighIter == rhs.neighIter);
    }
};

template <class GraphType, class NodeType, bool InEdges = false>
class NeighborRangeBase {
protected:
    const GraphType *G;
    NodeType source;

public:
    NeighborRangeBase(const GraphType &G, NodeType source) : G(&G), source(source) {
        assert(G.hasNode(source));
    }

    NeighborRangeBase() : G(nullptr) {}

    virtual typename GraphType::NeighborIterator begin() const = 0;
    virtual typename GraphType::NeighborIterator end() const = 0;
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_GENERAL_GRAPH_HPP_

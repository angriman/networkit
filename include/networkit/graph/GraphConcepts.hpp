/*  GraphConcepts.hpp
 *
 *  Created on: 05.07.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_GRAPH_GRAPH_CONCEPTS_HPP_
#define NETWORKIT_GRAPH_GRAPH_CONCEPTS_HPP_

#include <concepts>

#include <networkit/Globals.hpp>

namespace NetworKit {

template <typename T>
concept GraphNode = std::integral<T>;

template <typename T>
concept GraphEdgeWeight = std::integral<T> || std::floating_point<T>;

/**
 * Minimal read-only interface that algorithms, readers and writers templated on
 * the graph type rely on: the node and edge weight types, the sizes, and
 * whether the graph is directed or weighted.
 */
template <typename G>
concept GraphLike = requires(const G &g) {
    typename G::NodeT;
    typename G::EdgeWeightT;
    requires GraphNode<typename G::NodeT>;
    requires GraphEdgeWeight<typename G::EdgeWeightT>;
    { g.numberOfNodes() } -> std::convertible_to<count>;
    { g.numberOfEdges() } -> std::convertible_to<count>;
    { g.upperNodeIdBound() } -> std::convertible_to<typename G::NodeT>;
    { g.isDirected() } -> std::convertible_to<bool>;
    { g.isWeighted() } -> std::convertible_to<bool>;
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_GRAPH_CONCEPTS_HPP_

/*  GraphConcepts.hpp
 *
 *  Created on: 05.07.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_GRAPH_GRAPH_CONCEPTS_HPP_
#define NETWORKIT_GRAPH_GRAPH_CONCEPTS_HPP_

#include <concepts>

namespace NetworKit {

template <typename T>
concept GraphNode = std::integral<T>;

// Empty struct to represent an unweighted graph.
struct Unweighted {};

template <typename T>
concept GraphEdgeWeight = std::integral<T> || std::floating_point<T> || std::same_as<T, Unweighted>;

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_GRAPH_CONCEPTS_HPP_

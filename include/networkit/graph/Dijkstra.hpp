// networkit-format

#ifndef NETWORKIT_GRAPH_DIJKSTRA_HPP_
#define NETWORKIT_GRAPH_DIJKSTRA_HPP_

#include <array>
#include <limits>
#include <vector>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {
namespace Traversal {

/**
 * Iterate over nodes following Dijkstra's algorithm starting from the nodes in the given range.
 *
 * @param G The input graph.
 * @param first The first element of the range.
 * @param handle Takes a node as input parameter.
 */
template <class InputIt, typename L>
void DijkstraFrom(const Graph &G, InputIt first, InputIt last, L handle) {
    const node n = G.upperNodeIdBound();
    std::vector<edgeweight> distance(n, std::numeric_limits<edgeweight>::max());

    const auto lowerDistance = [&distance](node x, node y) noexcept->bool {
        return distance[x] < distance[y];
    };

    tlx::DAryAddressableIntHeap<node, 2, decltype(lowerDistance)> prioQ(lowerDistance);
    prioQ.reserve(n);
    for (; first != last; ++first) {
        distance[*first] = 0;
        prioQ.push(*first);
    }

    do {
        const node u = prioQ.extract_top();
        handle(u, distance[u]);

        node v;
        edgeweight weight;
        for (const auto neighWeight : G.weightNeighborRange(u)) {
            std::tie(v, weight) = neighWeight;
            const edgeweight newDist = distance[u] + weight;
            if (newDist < distance[v]) {
                distance[v] = newDist;
                prioQ.update(v);
            }
        }
    } while (!prioQ.empty());
}

/**
 * Iterate over nodes following Dijkstra's algorithm starting from the given source node.
 *
 * @param G The input graph.
 * @param source The source node.
 * @param handle Takes a node as input parameter.
 */
template <typename L>
void DijkstraFrom(const Graph &G, node source, L handle) {
    std::array<node, 1> startNodes{{source}};
    DijkstraFrom(G, startNodes.begin(), startNodes.end(), handle);
}

} // namespace Traversal
} // namespace NetworKit
#endif // NETWORKIT_GRAPH_DIJKSTRA_HPP_

#ifndef NETWORKIT_DISTANCE_BFS_IMPL_HPP_
#define NETWORKIT_DISTANCE_BFS_IMPL_HPP_

#include <limits>

#include <networkit/distance/SSSP.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit::BFSDetails {

/**
 * @ingroup distance
 * Runs a breadth-first search on the given graph from a given source node.
 */

template <bool Reverse = false>
class BFSImpl final : public SSSP {

    static constexpr auto infDist = std::numeric_limits<edgeweight>::max();

public:
    /*
     * @param G The graph
     * @param source The source node of the breadth-first search
     * @param storePaths Paths are reconstructable and the number of paths is
     * stored.
     * @param storeNodesSortedByDistance Store a vector of nodes ordered in
     * increasing distance from the source.
     * @param target The target node.
     */
    BFSImpl(const Graph &G, node source, bool storePaths = true,
            bool storeNodesSortedByDistance = false, node target = none)
        : SSSP(G, source, target, storePaths, storeNodesByDistance, target) {}

    void run() override;
};

} // namespace NetworKit::BFSDetails

#endif // NETWORKIT_DISTANCE_BFS_IMPL_HPP_

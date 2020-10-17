/*
 *  MultiSourceDijkstra.hpp
 *  Created on: 17.10.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#ifndef NETWORKIT_DISTANCE_MULTI_SOURCE_DIJKSTRA_HPP_
#define NETWORKIT_DISTANCE_MULTI_SOURCE_DIJKSTRA_HPP_

#include <networkit/distance/MSSP.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Multi-source Dijkstra search on a graph from a given set of source nodes to all the other nodes.
 * If one or more target nodes are specified, the Dijkstra terminates after the shortest-path
 * distance has been computed for all the target nodes.
 */
class MultiSourceDijkstra final : public MSSP {

    void runWithoutTargets();
    void runWithTargets();

public:
    /**
     * Creates an instance of MultiSourceDijkstra with the input graph @G.
     *
     * @param G The graph.
     */
    MultiSourceDijkstra(const Graph &G) : MSSP(G) {}

    /**
     * Creates an instance of MultiSourceDijkstra with the input graph @G and a given set of
     * sources.
     *
     * @param G The graph.
     * @param sourcesFirst,sourcesLast The range that contains the source nodes.
     */
    template <class SourcesInputIt>
    MultiSourceDijkstra(const Graph &G, SourcesInputIt sourcesFirst, SourcesInputIt sourcesLast)
        : MSSP(G, sourcesFirst, sourcesLast) {}

    /**
     * Creates an instance of MultiSourceDijkstra with the input graph @G and a given set of sources
     * and a given set of targets.
     *
     * @param G The graph.
     * @param sourcesFirst,sourcesLast The range that contains the source nodes.
     * @param targetsFirst,targetsLast The range that contains the targets nodes.
     */
    template <class SourcesInputIt, class TargetsInputIt>
    MultiSourceDijkstra(const Graph &G, SourcesInputIt sourcesFirst, SourcesInputIt sourcesLast,
                        TargetsInputIt targetsFirst, TargetsInputIt targetsLast)
        : MSSP(G, sourcesFirst, sourcesLast, targetsFirst, targetsLast) {}

    /**
     * Creates an instance of MultiSourceDijkstra with the input graph @G and a given set of sources
     * and a single target.
     *
     * @param G The graph.
     * @param sourcesFirst,sourcesLast The range that contains the source nodes.
     * @param target The target node.
     */
    template <class SourcesInputIt>
    MultiSourceDijkstra(const Graph &G, SourcesInputIt sourcesFirst, SourcesInputIt sourcesLast,
                        node target)
        : MSSP(G, sourcesFirst, sourcesLast, target) {}

    /*
     * Run the algorithm.
     */
    void run() override;

private:
    struct Less {
        Less(const std::vector<edgeweight> &dist) : dist(dist) {}
        bool operator()(node x, node y) const noexcept { return dist[x] < dist[y]; }

    private:
        const std::vector<edgeweight> &dist;
    };

    tlx::d_ary_addressable_int_heap<node, 2, Less> heap{dist};
};

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_MULTI_SOURCE_DIJKSTRA_HPP_

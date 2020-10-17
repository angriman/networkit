/*
 *  MSSP.hpp
 *  Created on: 17.10.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#ifndef NETWORKIT_DISTANCE_MSSP_HPP_
#define NETWORKIT_DISTANCE_MSSP_HPP_

#include <unordered_map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Abstract base class for multi-source shortest path algorithms.
 *
 * Computes the shortest-path distances from the source nodes to all the other nodes.
 * If one or more target nodes are specified, the algorithm terminates after the shortest-path
 * distance has been computed for all the target nodes.
 */
class MSSP : public Algorithm {

public:
    /**
     * Creates an instance of MSSP with the input graph @a G.
     *
     * @param G The graph.
     */
    MSSP(const Graph &G) : G(&G) {}

    /**
     * Creates an instance of MSSP with the input graph @a G and a given set of sources.
     *
     * @param G The graph.
     * @param sourcesFirst,sourcesLast The range that contains the source nodes.
     */
    template <class SourcesInputIt>
    MSSP(const Graph &G, SourcesInputIt sourcesFirst, SourcesInputIt sourcesLast) : MSSP(G) {
        setSources(sourcesFirst, sourcesLast);
    }

    /**
     * Creates an instance of MSSP with the input graph @a G and a given set of sources and a given
     * set of targets.
     *
     * @param G The graph.
     * @param sourcesFirst,sourcesLast The range that contains the source nodes.
     * @param targetsFirst,targetsLast The range that contains the target nodes.
     */
    template <class SourcesInputIt, class TargetsInputIt>
    MSSP(const Graph &G, SourcesInputIt sourcesFirst, SourcesInputIt sourcesLast,
         TargetsInputIt targetsFirst, TargetsInputIt targetsLast)
        : MSSP(G, sourcesFirst, sourcesLast) {
        setTargets(targetsFirst, targetsLast);
    }

    /**
     * Creates an instance of MSSP with the input graph @a G and a given set of sources and a single
     * target.
     *
     * @param G The graph.
     * @param sourcesFirst,sourcesLast The range that contains the source nodes.
     * @param target The target node.
     */
    template <class SourcesInputIt>
    MSSP(const Graph &G, SourcesInputIt sourcesFirst, SourcesInputIt sourcesLast, node target)
        : MSSP(G, sourcesFirst, sourcesLast) {
        setTarget(target);
    }

    ~MSSP() override = default;

    /**
     * Run the algorithm.
     */
    void run() override = 0;

    /**
     * Returns the minimum distances from the source nodes to all the other nodes. If one or more
     * targets have been specified, only the distances to the targets are guaranteed to be correct.
     *
     * @return The distances from the source nodes to all the other nodes.
     */
    const std::vector<edgeweight> &getDistances() const;

    /**
     * Returns the distance from the source nodes to the node @a t. If one or more targets have been
     * specified and @a t is not among them, the output distance is not guaranteed to be correct.
     *
     * @param t A node.
     * @return The minimum distance from the source nodes to the node @a t.
     */
    edgeweight distance(node t) const;

    /**
     * Returns a vector with all the nodes of the graph sorted by increasing distance from the
     * source nodes.
     *
     * @return Vector of nodes sorted by increasing distance from the source nodes.
     */
    std::vector<node> getNodesSortedByDistance() const;

    /**
     * Returns a vector with the target nodes sorted by increasing distance from the source nodes.
     *
     * @return Vector of target nodes sorted by increasing distance from the source nodes.
     */
    std::vector<node> getTargetNodesSortedByDistance() const;

    /**
     * Returns a <node, edgeweight> map with the minimum distances from the source nodes to the
     * target nodes.
     *
     * @return A <node, edgeweight> map with keys the target nodes and values their minimum
     * distances from the source nodes.
     */
    std::unordered_map<node, edgeweight> getDistancesToTargets() const;

    /**
     * Returns the number of nodes reached from the source nodes during the exploration of the
     * graph.
     *
     * @return The number of reached nodes from the source nodes during the exploration of the
     * graph.
     */
    count getReachableNodes() const;

    /**
     * Change the source nodes.
     *
     * @param sourcesFirst,sourcesLast The range that contains the new set of source nodes.
     */
    template <class SourcesInputIt>
    void setSources(SourcesInputIt sourcesFirst, SourcesInputIt sourcesLast) {
        sources.clear();
        sources.insert(sources.begin(), sourcesFirst, sourcesLast);
    }

    /**
     * Change the source nodes.
     *
     * @param source The new source node.
     */
    void setSource(node source) { sources = {source}; }

    /**
     * Change the target nodes.
     *
     * @param targetsFirst,targetsLast The range that contains the new set of target nodes.
     */
    template <class TargetsInputIt>
    void setTargets(TargetsInputIt targetsFirst, TargetsInputIt targetsLast) {
        targets.clear();
        targets.insert(targets.begin(), targetsFirst, targetsLast);
    }

    /**
     * Change the target nodes.
     *
     * @param target The new target node.
     */
    void setTarget(node target) { targets = {target}; }

    /**
     * Clear the target nodes.
     */
    void clearTargets() { targets.clear(); }

protected:
    const Graph *G;
    std::vector<node> sources, targets;
    count reachedNodes;
    std::vector<edgeweight> dist;

    static constexpr edgeweight infdist = std::numeric_limits<edgeweight>::max();
};

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_MSSP_HPP_

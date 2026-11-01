// networkit-format
#ifndef NETWORKIT_GRAPH_MINIMUM_ARBORESCENCE_HPP_
#define NETWORKIT_GRAPH_MINIMUM_ARBORESCENCE_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 * Computes a minimum weight arborescence of a directed graph, i.e., the lightest set of edges that
 * contains exactly one path from a given root to every node. This is the directed counterpart of a
 * minimum spanning tree, which greedy algorithms cannot compute. The Chu-Liu/Edmonds algorithm
 * runs in O(nm) time and space.
 *
 * Unweighted graphs are treated as having unit weights, and negative weights are allowed.
 * Self-loops and edges into the root are ignored.
 */
class MinimumArborescence final : public Algorithm {

public:
    /**
     * Constructs the algorithm for the Graph @a G.
     *
     * @param[in] G A directed graph.
     * @param[in] root The node all paths of the arborescence start from.
     * @throws std::runtime_error if @a G is undirected or does not contain @a root.
     */
    MinimumArborescence(const Graph &G, node root);

    /**
     * Computes the minimum arborescence of the graph passed in the constructor.
     */
    void run() override;

    /**
     * Returns whether an arborescence exists, i.e., whether the root reaches all nodes.
     *
     * @return True if the graph has an arborescence with the given root.
     */
    bool hasArborescence() const {
        assureFinished();
        return exists;
    }

    /**
     * Returns the minimum arborescence, a directed and weighted graph with the same node ids as the
     * input graph. It has no edges if there is no arborescence.
     *
     * @return The minimum arborescence.
     */
    Graph getArborescence() const {
        assureFinished();
        return arborescence;
    }

    /**
     * Returns the total weight of the minimum arborescence.
     *
     * @return The total weight of the edges of the arborescence.
     */
    edgeweight getTotalWeight() const {
        assureFinished();
        return totalWeight;
    }

private:
    const Graph *G;
    node root;
    bool exists = false;
    edgeweight totalWeight = 0;
    Graph arborescence;
};

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_MINIMUM_ARBORESCENCE_HPP_

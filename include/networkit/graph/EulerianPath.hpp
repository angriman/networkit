// networkit-format
#ifndef NETWORKIT_GRAPH_EULERIAN_PATH_HPP_
#define NETWORKIT_GRAPH_EULERIAN_PATH_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 * Computes an Eulerian path of a graph, i.e., a walk that uses every edge exactly once, with
 * Hierholzer's algorithm in O(n + m). If such a walk can start and end at the same node, the result
 * is a circuit. Edge weights are ignored; parallel edges and self-loops are each used once.
 * Nodes without edges are ignored, so they do not make the graph disconnected.
 */
class EulerianPath final : public Algorithm {

public:
    /**
     * Constructs the algorithm for the Graph @a G.
     *
     * @param[in] G Directed or undirected graph.
     */
    explicit EulerianPath(const Graph &G);

    /**
     * Computes an Eulerian path of the graph passed in the constructor, if there is one.
     */
    void run() override;

    /**
     * Returns whether the graph has an Eulerian path. This is the case for graphs without edges.
     *
     * @return True if the graph has an Eulerian path.
     */
    bool hasEulerianPath() const {
        assureFinished();
        return exists;
    }

    /**
     * Returns the nodes along the Eulerian path, so that consecutive nodes are connected by an edge
     * and the path has m + 1 nodes. It is empty if there is no Eulerian path or the graph has no
     * edges.
     *
     * @return The nodes along the Eulerian path.
     */
    const std::vector<node> &getPath() const {
        assureFinished();
        return path;
    }

private:
    const Graph *G;
    bool exists = false;
    std::vector<node> path;
};

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_EULERIAN_PATH_HPP_

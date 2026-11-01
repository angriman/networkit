// networkit-format
#ifndef NETWORKIT_FLOW_STOER_WAGNER_HPP_
#define NETWORKIT_FLOW_STOER_WAGNER_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup flow
 * Computes a global minimum cut of an undirected graph, i.e., a partition of its nodes into two
 * non-empty sets such that the total weight of the edges between them is minimal, with the
 * algorithm from "A simple min-cut algorithm" by M. Stoer and F. Wagner (1997). Unlike a maximum
 * flow, this requires neither a source nor a target.
 *
 * Unweighted graphs are treated as having unit weights. Self-loops are ignored and the weights of
 * parallel edges are summed up. The running time is in O(nm log n).
 */
class StoerWagner final : public Algorithm {

public:
    /**
     * Constructs the algorithm for the Graph @a G.
     *
     * @param G An undirected graph with at least two nodes and non-negative edge weights.
     * @throws std::runtime_error if @a G is directed, has fewer than two nodes, or has an edge with
     * negative weight.
     */
    explicit StoerWagner(const Graph &G);

    /**
     * Computes a global minimum cut of the graph passed in the constructor.
     */
    void run() override;

    /**
     * Returns the total weight of the edges of the minimum cut.
     *
     * @return The weight of the minimum cut.
     */
    edgeweight getCutValue() const {
        assureFinished();
        return cutValue;
    }

    /**
     * Returns the nodes on one side of the minimum cut, in ascending order. The other side consists
     * of all the remaining nodes of the graph.
     *
     * @return The nodes on one side of the minimum cut.
     */
    const std::vector<node> &getCutSide() const {
        assureFinished();
        return cutSide;
    }

private:
    const Graph *G;
    edgeweight cutValue = 0;
    std::vector<node> cutSide;
};

} // namespace NetworKit

#endif // NETWORKIT_FLOW_STOER_WAGNER_HPP_

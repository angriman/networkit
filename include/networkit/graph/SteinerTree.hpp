// networkit-format
#ifndef NETWORKIT_GRAPH_STEINER_TREE_HPP_
#define NETWORKIT_GRAPH_STEINER_TREE_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 * Computes a light tree that connects a given set of terminal nodes, possibly through other nodes.
 * Finding the lightest one is NP-hard. This is the approximation algorithm from "A faster
 * approximation algorithm for the Steiner problem in graphs" by K. Mehlhorn (1988): the weight of
 * the result is at most 2(1 - 1/l) times the optimal one, for l terminals. It runs in
 * O(m log n) time.
 *
 * Unweighted graphs are treated as having unit weights.
 */
class SteinerTree final : public Algorithm {

public:
    /**
     * Constructs the algorithm for the Graph @a G.
     *
     * @param[in] G An undirected graph with non-negative edge weights.
     * @param[in] terminals The nodes to connect. Duplicates are ignored.
     * @throws std::runtime_error if @a G is directed or has negative edge weights, or if one of
     * the terminals is not in @a G.
     */
    SteinerTree(const Graph &G, std::vector<node> terminals);

    /**
     * Computes the tree for the graph and terminals passed in the constructor.
     */
    void run() override;

    /**
     * Returns whether all terminals are in the same connected component.
     *
     * @return True if a tree connecting all terminals exists.
     */
    bool hasTree() const {
        assureFinished();
        return connected;
    }

    /**
     * Returns the tree, a weighted graph with the same node ids as the input graph. Nodes outside
     * the tree have no edges. It has no edges if the terminals are not connected.
     *
     * @return The tree connecting the terminals.
     */
    Graph getTree() const {
        assureFinished();
        return tree;
    }

    /**
     * Returns the total weight of the tree.
     *
     * @return The total weight of the edges of the tree.
     */
    edgeweight getTotalWeight() const {
        assureFinished();
        return totalWeight;
    }

private:
    const Graph *G;
    std::vector<node> terminals;
    bool connected = false;
    edgeweight totalWeight = 0;
    Graph tree;
};

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_STEINER_TREE_HPP_

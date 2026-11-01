// networkit-format
#ifndef NETWORKIT_FLOW_GOMORY_HU_TREE_HPP_
#define NETWORKIT_FLOW_GOMORY_HU_TREE_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup flow
 * Computes a Gomory-Hu tree of an undirected graph with the algorithm from "Very simple methods for
 * all pairs network flow analysis" by D. Gusfield (1990). The tree has the same nodes as the graph,
 * and the minimum cut between two nodes equals the lightest edge on their tree path. Removing that
 * edge splits the nodes into the two sides of such a cut. This answers all minimum cut queries with
 * n - 1 maximum flow computations instead of one per pair.
 *
 * Unweighted graphs are treated as having unit weights. Self-loops are ignored and the weights of
 * parallel edges are summed up.
 */
class GomoryHuTree final : public Algorithm {

public:
    /**
     * Constructs the algorithm for the Graph @a G.
     *
     * @param[in] G An undirected graph with non-negative edge weights.
     * @throws std::runtime_error if @a G is directed or has an edge with negative weight.
     */
    explicit GomoryHuTree(const Graph &G);

    /**
     * Computes the Gomory-Hu tree of the graph passed in the constructor.
     */
    void run() override;

    /**
     * Returns the Gomory-Hu tree, a weighted graph with the same node ids as the input graph.
     * Nodes in different components of the input are connected by edges of weight 0.
     *
     * @return The Gomory-Hu tree.
     */
    Graph getTree() const {
        assureFinished();
        return tree;
    }

    /**
     * Returns the value of the minimum cut between two nodes in time linear in their tree distance.
     *
     * @param[in] u A node of the graph.
     * @param[in] v Another node of the graph.
     * @return The weight of a minimum cut separating @a u and @a v.
     * @throws std::runtime_error if @a u and @a v are equal or not in the graph.
     */
    edgeweight getMinCutValue(node u, node v) const;

private:
    const Graph *G;
    Graph tree;
    std::vector<node> parent;
    std::vector<edgeweight> parentWeight;
    std::vector<count> depth;
};

} // namespace NetworKit

#endif // NETWORKIT_FLOW_GOMORY_HU_TREE_HPP_

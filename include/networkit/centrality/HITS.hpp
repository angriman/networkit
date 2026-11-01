// networkit-format
#ifndef NETWORKIT_CENTRALITY_HITS_HPP_
#define NETWORKIT_CENTRALITY_HITS_HPP_

#include <vector>

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * Computes the hub and authority scores of the HITS algorithm from "Authoritative sources in a
 * hyperlinked environment" by J. Kleinberg (1999). A node is a good authority if good hubs point
 * to it, and a good hub if it points to good authorities. Edge weights multiply the contribution
 * of an edge.
 *
 * Both score vectors have Euclidean norm 1, unless the graph has no edges, in which case all scores
 * are 0. On connected, undirected, non-bipartite graphs hub and authority scores coincide; on
 * bipartite graphs they may differ.
 */
class HITS final : public Centrality {

public:
    /**
     * Constructs the HITS class for the Graph @a G.
     *
     * @param[in] G Graph to be processed.
     * @param[in] tol Error tolerance: the algorithm stops when the Euclidean distance between the
     * scores of two consecutive iterations is at most @a tol. Must be positive.
     * @param[in] maxIterations Maximum number of iterations. Must be positive.
     */
    HITS(const Graph &G, double tol = 1e-8, count maxIterations = 1000);

    /**
     * Computes the hub and authority scores of the graph passed in the constructor.
     */
    void run() override;

    /**
     * Returns the maximum possible score, reached by a node that holds the whole Euclidean norm.
     */
    double maximum() override { return 1.; }

    /**
     * Returns the authority score of each node, same as scores().
     *
     * @return Authority scores indexed by node id.
     */
    const std::vector<double> &authorityScores() const { return scores(); }

    /**
     * Returns the hub score of each node.
     *
     * @return Hub scores indexed by node id.
     */
    const std::vector<double> &hubScores() const {
        assureFinished();
        return hubData;
    }

    /**
     * Returns the number of iterations performed by the algorithm.
     *
     * @return Number of iterations performed by the algorithm.
     */
    count numberOfIterations() const {
        assureFinished();
        return iterations;
    }

private:
    double tol;
    count maxIterations;
    count iterations = 0;
    std::vector<double> hubData;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_HITS_HPP_

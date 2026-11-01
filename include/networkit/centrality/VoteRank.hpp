// networkit-format
#ifndef NETWORKIT_CENTRALITY_VOTE_RANK_HPP_
#define NETWORKIT_CENTRALITY_VOTE_RANK_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * Identifies a set of influential spreaders with the VoteRank algorithm from "Identifying a set of
 * influential spreaders in complex networks" by J. Zhang et al. (2016). Nodes vote for their
 * neighbors and the node with the most votes is elected in each round. The elected node stops
 * voting and its neighbors lose part of their voting ability, which keeps the spreaders apart.
 *
 * In directed graphs a node receives the votes of its out-neighbors, and an elected node weakens
 * them. Edge weights are ignored. The running time is in O(k(n + m)) for k spreaders.
 */
class VoteRank final : public Algorithm {

public:
    /**
     * Constructs the VoteRank class for the Graph @a G.
     *
     * @param[in] G Graph to be processed.
     * @param[in] k Maximum number of spreaders. The algorithm stops earlier if no node gets votes.
     * By default it continues until then.
     */
    explicit VoteRank(const Graph &G, count k = none);

    /**
     * Elects the spreaders of the graph passed in the constructor.
     */
    void run() override;

    /**
     * Returns the elected spreaders, in order of election.
     *
     * @return The spreaders, the most influential one first.
     */
    const std::vector<node> &getSpreaders() const {
        assureFinished();
        return spreaders;
    }

private:
    const Graph *G;
    count k;
    std::vector<node> spreaders;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_VOTE_RANK_HPP_

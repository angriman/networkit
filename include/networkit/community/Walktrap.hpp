// networkit-format
#ifndef NETWORKIT_COMMUNITY_WALKTRAP_HPP_
#define NETWORKIT_COMMUNITY_WALKTRAP_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Walktrap community detection from "Computing communities in large networks using random walks"
 * by P. Pons and M. Latapy (2005). Short random walks tend to stay in the community they start
 * from, so nodes that reach the rest of the graph with similar probabilities are likely to belong
 * together. The algorithm starts from singletons, repeatedly merges the two adjacent communities
 * with the most similar walks, and returns the partition with the highest modularity seen.
 *
 * Connected components are never merged. Walks are stored as dense vectors, so the algorithm needs
 * O(n^2) memory and O(mn(t + H)) time for walks of length t and a dendrogram of height H, which is
 * O(mn^2) in the worst case. This limits it to graphs of up to some ten thousand nodes.
 *
 * Each merge uses the exact costs of all candidate pairs. The authors' implementation, which is
 * also the one in igraph, estimates most of them to save time, so its results can differ.
 */
class Walktrap final : public CommunityDetectionAlgorithm {

public:
    /**
     * Constructs the algorithm for the Graph @a G.
     *
     * @param[in] G An undirected graph with non-negative edge weights.
     * @param[in] walkLength The length of the random walks. Must be positive.
     * @throws std::runtime_error if @a G is directed or has negative edge weights, or if
     * @a walkLength is zero.
     */
    explicit Walktrap(const Graph &G, count walkLength = 4);

    /**
     * Detects the communities of the graph passed in the constructor.
     */
    void run() override;

private:
    count walkLength;
};

} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_WALKTRAP_HPP_

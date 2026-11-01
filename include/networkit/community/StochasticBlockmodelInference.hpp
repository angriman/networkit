// networkit-format
#ifndef NETWORKIT_COMMUNITY_STOCHASTIC_BLOCKMODEL_INFERENCE_HPP_
#define NETWORKIT_COMMUNITY_STOCHASTIC_BLOCKMODEL_INFERENCE_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Partitions the nodes of an undirected graph into k blocks by fitting a stochastic blockmodel,
 * using the profile likelihood of Karrer and Newman, "Stochastic blockmodels and community
 * structure in networks", Physical Review E 83, 016107 (2011).
 *
 * Unlike modularity-based methods, the model does not assume assortative communities: blocks can
 * also be, for example, the two sides of a bipartite graph. The degree-corrected variant keeps
 * nodes of very different degrees in the same block, which the plain variant tends to separate.
 *
 * Edge weights are read as the means of Poisson-distributed edge multiplicities, so they must not
 * be negative. The search maximizes the likelihood with local moves from random starting points,
 * so the result is a local optimum. Restarts run in parallel, and the result is deterministic for
 * a fixed seed only when running on a single thread.
 *
 * With s the average number of blocks adjacent to a node, a sweep takes O(m + n k s) time. The
 * algorithm needs O(n + m + k^2) memory per thread.
 */
class StochasticBlockmodelInference final : public CommunityDetectionAlgorithm {
public:
    /**
     * @param G Undirected graph, with non-negative edge weights.
     * @param k Number of blocks, between 1 and the number of nodes. Blocks may end up empty.
     * @param degreeCorrected Whether to use the degree-corrected variant of the model.
     * @param restarts Number of random starting points, at least 1.
     * @param maxSweeps Maximum number of passes over the nodes per starting point, at least 1.
     * @throws std::runtime_error if G is directed.
     * @throws std::invalid_argument if G has negative edge weights, or if k, restarts or maxSweeps
     * are out of range.
     */
    StochasticBlockmodelInference(const Graph &G, count k, bool degreeCorrected = true,
                                  count restarts = 30, count maxSweeps = 100);

    void run() override;

    /**
     * Returns the profile log-likelihood of the partition, up to a constant that does not depend
     * on the partition.
     */
    double getLogLikelihood() const;

private:
    count k;
    bool degreeCorrected;
    count restarts;
    count maxSweeps;
    double logLikelihood = 0.;
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_STOCHASTIC_BLOCKMODEL_INFERENCE_HPP_

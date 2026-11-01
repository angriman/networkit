// networkit-format
#ifndef NETWORKIT_MATCHING_BLOSSOM_MATCHER_HPP_
#define NETWORKIT_MATCHING_BLOSSOM_MATCHER_HPP_

#include <networkit/matching/Matcher.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * Computes a maximum weight matching of a general graph with the primal-dual blossom method from
 * "Paths, trees, and flowers" by J. Edmonds (1965). Unlike the other matchers for general graphs,
 * which are approximations, the result is optimal. The running time is in O(n^2 (n + m log n)) in
 * the worst case, and much lower in practice.
 *
 * Edges without a positive weight are never matched, and the matching does not have to be perfect
 * or of maximum cardinality. On unweighted graphs, where all weights are one, it is a maximum
 * cardinality matching. Self-loops are ignored, and only the heaviest of parallel edges is used.
 * Weights are floating-point numbers, so results for weights that are not integers are optimal up
 * to rounding errors.
 */
class BlossomMatcher final : public Matcher {

public:
    /**
     * Constructs the matcher for the Graph @a G.
     *
     * @param[in] G An undirected graph.
     * @throws std::invalid_argument if @a G is directed.
     */
    explicit BlossomMatcher(const Graph &G);

    /**
     * Computes a maximum weight matching of the graph passed in the constructor.
     */
    void run() override;
};

} /* namespace NetworKit */

#endif // NETWORKIT_MATCHING_BLOSSOM_MATCHER_HPP_

// networkit-format
#ifndef NETWORKIT_MATCHING_HOPCROFT_KARP_MATCHER_HPP_
#define NETWORKIT_MATCHING_HOPCROFT_KARP_MATCHER_HPP_

#include <networkit/matching/Matcher.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * Computes a maximum cardinality matching of a bipartite graph with the algorithm from "An n^(5/2)
 * algorithm for maximum matchings in bipartite graphs" by J. Hopcroft and R. Karp (1973). The
 * running time is in O(m sqrt(n)). Edge weights are ignored.
 */
class HopcroftKarpMatcher final : public Matcher {

public:
    /**
     * Constructs the matcher for the Graph @a G.
     *
     * @param[in] G An undirected, bipartite graph. The two sides are detected automatically.
     * @throws std::invalid_argument if @a G is directed or not bipartite, which includes graphs
     * with self-loops.
     */
    explicit HopcroftKarpMatcher(const Graph &G);

    /**
     * Computes a maximum cardinality matching of the graph passed in the constructor.
     */
    void run() override;
};

} /* namespace NetworKit */

#endif // NETWORKIT_MATCHING_HOPCROFT_KARP_MATCHER_HPP_

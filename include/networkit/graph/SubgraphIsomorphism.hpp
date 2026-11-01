// networkit-format
#ifndef NETWORKIT_GRAPH_SUBGRAPH_ISOMORPHISM_HPP_
#define NETWORKIT_GRAPH_SUBGRAPH_ISOMORPHISM_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 * Finds the occurrences of a pattern graph in a target graph with a backtracking search in the
 * style of "A (sub)graph isomorphism algorithm for matching large graphs" by L. Cordella et al.
 * (2004). A match assigns a different target node to each pattern node, such that each pattern edge
 * is mapped to a target edge. An induced match additionally forbids target edges between the
 * matched nodes that are not in the pattern.
 *
 * Symmetric patterns are matched once per automorphism: a triangle has six matches on each triangle
 * of the target. Edge weights and the multiplicity of parallel edges are ignored. The worst-case
 * running time is exponential in the size of the pattern.
 */
class SubgraphIsomorphism final : public Algorithm {

public:
    /**
     * Constructs the search for the pattern @a P in the Graph @a G.
     *
     * @param[in] G The target graph.
     * @param[in] P The pattern to look for. Must have at least one node.
     * @param[in] induced If true, matched nodes must not have additional edges among them.
     * @param[in] maxMatches Stop after this number of matches. By default, all matches are stored,
     * of which there can be many more than target nodes.
     * @throws std::invalid_argument if @a P has no nodes, @a maxMatches is 0, or only one of the
     * graphs is directed.
     */
    SubgraphIsomorphism(const Graph &G, const Graph &P, bool induced = true,
                        count maxMatches = none);

    /**
     * Searches for the matches.
     */
    void run() override;

    /**
     * Returns the matches found. Each one is indexed by pattern node id and contains the matched
     * target node, or none for ids that are not in the pattern.
     *
     * @return The matches found.
     */
    const std::vector<std::vector<node>> &getMatches() const {
        assureFinished();
        return matches;
    }

    /**
     * Returns the number of matches found.
     *
     * @return The number of matches found.
     */
    count numberOfMatches() const {
        assureFinished();
        return matches.size();
    }

private:
    const Graph *G;
    const Graph *P;
    bool induced;
    count maxMatches;
    std::vector<std::vector<node>> matches;
};

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_SUBGRAPH_ISOMORPHISM_HPP_

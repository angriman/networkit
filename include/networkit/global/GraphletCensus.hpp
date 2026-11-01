// networkit-format
#ifndef NETWORKIT_GLOBAL_GRAPHLET_CENSUS_HPP_
#define NETWORKIT_GLOBAL_GRAPHLET_CENSUS_HPP_

#include <array>
#include <cstdint>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup global
 * Counts the connected undirected graphlets on three and four nodes, i.e., the sets of nodes that
 * induce each of these subgraphs. A triangle with a pendant edge is a paw, a clique that misses an
 * edge is a diamond. For the directed motifs on three nodes, see TriadCensus.
 *
 * The counts are derived from subgraph counts with closed formulas rather than by visiting each
 * occurrence, which only cliques require on top of the triangles. The algorithm runs in O(m d^2)
 * time for a maximum degree d. Edge weights, the multiplicity of parallel edges and self-loops are
 * ignored.
 */
class GraphletCensus final : public Algorithm {

public:
    enum class Graphlet : uint8_t {
        Wedge,    // Path on 3 nodes
        Triangle, // Cycle on 3 nodes
        Path,     // Path on 4 nodes
        Star,     // A node with 3 neighbors
        Cycle,    // Cycle on 4 nodes
        Paw,      // Triangle with a pendant edge
        Diamond,  // Clique on 4 nodes without an edge
        Clique    // Clique on 4 nodes
    };

    static constexpr count numberOfGraphlets = 8;

    /**
     * Constructs the algorithm for the Graph @a G.
     *
     * @param[in] G An undirected graph.
     * @throws std::runtime_error if @a G is directed.
     */
    explicit GraphletCensus(const Graph &G);

    /**
     * Computes the census of the graph passed in the constructor.
     */
    void run() override;

    /**
     * Returns the number of sets of nodes that induce the given graphlet.
     *
     * @param[in] graphlet The graphlet to look up.
     * @return The number of its occurrences.
     */
    count getCount(Graphlet graphlet) const {
        assureFinished();
        return counts[static_cast<index>(graphlet)];
    }

    /**
     * Returns the number of occurrences of all graphlets, indexed in the order of Graphlet.
     *
     * @return The census of the graph.
     */
    std::vector<count> getCounts() const {
        assureFinished();
        return {counts.begin(), counts.end()};
    }

private:
    const Graph *G;
    std::array<count, numberOfGraphlets> counts{};
};

} /* namespace NetworKit */

#endif // NETWORKIT_GLOBAL_GRAPHLET_CENSUS_HPP_

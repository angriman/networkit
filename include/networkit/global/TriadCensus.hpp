// networkit-format
#ifndef NETWORKIT_GLOBAL_TRIAD_CENSUS_HPP_
#define NETWORKIT_GLOBAL_TRIAD_CENSUS_HPP_

#include <array>
#include <cstdint>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup global
 * Counts the directed motifs on three nodes: every triple of nodes forms one of 16 triads, and the
 * census tells how many triples form each of them. Comparing it with the census of a random graph
 * reveals over-represented motifs, such as feed-forward loops. This is the algorithm from "A
 * subquadratic triad census algorithm for large sparse networks with small maximum degree" by
 * V. Batagelj and A. Mrvar (2001), which only visits connected triples and runs in
 * O(m d log d) time for a maximum degree d.
 *
 * Edge weights, the multiplicity of parallel edges and self-loops are ignored.
 */
class TriadCensus final : public Algorithm {

public:
    /**
     * The 16 triads in the standard MAN notation: the digits are the numbers of mutual, asymmetric
     * and null dyads, and the letter tells apart triads where these are the same (Down, Up,
     * Cyclic, Transitive).
     */
    enum class Triad : uint8_t {
        T003,  // A  B  C
        T012,  // A->B  C
        T102,  // A<->B  C
        T021D, // A<-B->C
        T021U, // A->B<-C
        T021C, // A->B->C
        T111D, // A<->B<-C
        T111U, // A<->B->C
        T030T, // A->B<-C, A->C
        T030C, // A->B->C->A
        T201,  // A<->B<->C
        T120D, // A<-B->C, A<->C
        T120U, // A->B<-C, A<->C
        T120C, // A->B->C, A<->C
        T210,  // A->B<->C, A<->C
        T300   // A<->B<->C, A<->C
    };

    static constexpr count numberOfTriads = 16;

    /**
     * Constructs the algorithm for the Graph @a G.
     *
     * @param[in] G A directed graph.
     * @throws std::runtime_error if @a G is undirected.
     */
    explicit TriadCensus(const Graph &G);

    /**
     * Computes the census of the graph passed in the constructor.
     */
    void run() override;

    /**
     * Returns the number of triples of nodes that form the given triad.
     *
     * @param[in] triad The triad to look up.
     * @return The number of its occurrences.
     */
    count getCount(Triad triad) const {
        assureFinished();
        return counts[static_cast<index>(triad)];
    }

    /**
     * Returns the number of occurrences of all triads, indexed in the order of Triad.
     *
     * @return The census of the graph.
     */
    std::vector<count> getCounts() const {
        assureFinished();
        return {counts.begin(), counts.end()};
    }

private:
    const Graph *G;
    std::array<count, numberOfTriads> counts{};
};

} /* namespace NetworKit */

#endif // NETWORKIT_GLOBAL_TRIAD_CENSUS_HPP_

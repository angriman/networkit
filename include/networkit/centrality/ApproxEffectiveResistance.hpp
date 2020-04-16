/*
 * ApproxEffectiveResistance.hpp
 *
 *  Created on: 17.10.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#ifndef NETWORKIT_CENTRALITY_APPROX_EFFECTIVE_RESISTANCE_HPP_
#define NETWORKIT_CENTRALITY_APPROX_EFFECTIVE_RESISTANCE_HPP_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/BiconnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

#define PCG32_MULT 0x5851f42d4c957f2dULL

namespace NetworKit {

using small_node = uint32_t;
static constexpr small_node inf = std::numeric_limits<small_node>::max();

// See https://github.com/wjakob/pcg32
struct pcg32 {
    /// Initialize the pseudorandom number generator with the \ref seed() function
    pcg32(uint64_t initstate, uint64_t initseq = 1u) { seed(initstate, initseq); }

    void seed(uint64_t initstate, uint64_t initseq = 1) {
        state = 0U;
        inc = (initseq << 1u) | 1u;
        nextUInt();
        state += initstate;
        nextUInt();
    }

    /// Generate a uniformly distributed unsigned 32-bit random number
    uint32_t nextUInt() {
        uint64_t oldstate = state;
        state = oldstate * PCG32_MULT + inc;
        uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
        uint32_t rot = (uint32_t)(oldstate >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
    }

    /// Generate a uniformly distributed number, r, where 0 <= r < bound
    uint32_t nextUInt(uint32_t bound) {

        uint32_t threshold = (~bound + 1u) % bound;

        for (;;) {
            uint32_t r = nextUInt();
            if (r >= threshold)
                return r % bound;
        }
    }

    uint64_t state; // RNG state.  All values are possible.
    uint64_t inc;   // Controls which RNG sequence (stream) is selected. Must *always* be odd.
};

// NOTE: Does not support multiple calls of run()
class ApproxEffectiveResistance final : public Algorithm {

public:
    ApproxEffectiveResistance(const Graph &G, double epsilon = 0.1);

    ~ApproxEffectiveResistance() = default;

    void run() override;

    small_node getRoot() const noexcept { return root; }
    small_node getRootEcc() const noexcept { return rootEcc; }

    std::vector<int> getNonNormalizedData() const {
        std::vector<int> aggregated(G.upperNodeIdBound(), 0);
        G.parallelForNodes([&](const small_node u) {
            for (const auto &threadScores : approxEffResistanceGlobal) {
                aggregated[u] += threadScores[u];
            }
        });
        return aggregated;
    }

    // Aggregates per-thread scores into a single vector and returns it
    std::vector<double> getApproxEffectiveResistances() const {
        std::vector<double> result(G.upperNodeIdBound());
        G.parallelForNodes([&](const small_node u) {
            for (const auto &threadScores : approxEffResistanceGlobal) {
                result[u] += static_cast<double>(threadScores[u]);
            }
            result[u] /= static_cast<double>(numberOfUSTs);
        });

        return result;
    }

    count getRootEccentricity() const { return rootEcc; }
    count computeNumberOfUSTs() const {
        return rootEcc * rootEcc
               * static_cast<count>(
                   std::ceil(std::log(2.0 * static_cast<double>(G.numberOfEdges()) / delta)
                             / (2.0 * epsilon * epsilon)));
    }

    count numberOfUSTs = 0;
    double getEpsilon() const { return epsilon; }
    double timeToSample, timeDFS, timeToAggregate;

    void init();

private:
    // Input parameters
    const Graph &G;
    const double epsilon, delta;
    count sampledUSTs = 0;
    small_node root;
    uint32_t rootEcc;
    bool didInit = false;
    std::vector<count> samplingTime, dfsTime, aggregationTime;

    static constexpr int sweeps = 10;
    enum class NodeStatus : unsigned char {
        NOT_IN_COMPONENT,
        IN_SPANNING_TREE,
        NOT_VISITED,
        VISITED
    };

    // Used to mark the status of each node, one vector per thread
    std::vector<std::vector<NodeStatus>> statusGlobal;

    std::unique_ptr<BiconnectedComponents> bcc;

    // Nodes in each biconnected components sorted by their degree.
    std::vector<std::vector<small_node>> sequences;

    // Pointers to the parent of the UST, one vector per thread
    std::vector<std::vector<small_node>> parentGlobal;

    // Index of the parent component of the current component (after the topological order has been
    // determined)
    std::vector<index> biParent;
    // Node within the bionnected component that points to the node in the parent component
    std::vector<small_node> biAnchor;

    // Topological order of the biconencted components
    std::vector<index> topOrder;

    // Non-normalized approx effective resistance, one per thread.
    std::vector<std::vector<int>> approxEffResistanceGlobal;

    // Timestamps for DFS
    std::vector<std::vector<uint32_t>> tVisitGlobal, tFinishGlobal;

    // Random number generators
    std::vector<pcg32> generators;

    // Parent pointers of the bfs tree
    std::vector<small_node> bfsParent;

    // Nodes sequences: Wilson's algorithm runs faster if we start the random walks following a
    // specific sequence of nodes. In this function we compute those sequences.
    void computeNodeSequence();

    // Adjacency list for trees: additional data structure to speed-up the DFS
    std::vector<std::vector<small_node>> ustChildPtrGlobal;
    std::vector<std::vector<small_node>> ustSiblingPtrGlobal;

    void computeBFSTree();
    void sampleUST();
    void dfsUST();
    void aggregateUST();

    small_node approxMinEccNode();

    // Debugging methods
    void checkBFSTree() const;
    void checkUST() const;
    void checkTwoNodesSequence(const std::vector<small_node> &sequence) const;
    void checkTimeStamps() const;
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_APPROX_EFFECTIVE_RESISTANCE_HPP_

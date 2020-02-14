/*
 * ApproxEffectiveResistance.hpp
 *
 *  Created on: 17.10.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/BiconnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/numerics/ConjugateGradient.hpp>
#include <networkit/numerics/Preconditioner/IdentityPreconditioner.hpp>

#define PCG32_MULT 0x5851f42d4c957f2dULL

namespace NetworKit {

enum RootStrategy { MaxDegree, Random, MinApproxEcc };

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
    ApproxEffectiveResistance(const Graph &G, double epsilon = 0.1, double tolerance = 1e-6);

    ~ApproxEffectiveResistance() = default;

    void run() override;

    node getRoot() const noexcept { return root; }
    node getRootEcc() const noexcept { return rootEcc; }

    std::vector<int> getNonNormalizedData() const {
        std::vector<int> aggregated(G.upperNodeIdBound(), 0);
        G.parallelForNodes([&](const node u) {
            for (const auto &threadScores : approxEffResistanceGlobal) {
                aggregated[u] += threadScores[u];
            }
        });
        return aggregated;
    }

    // Aggregates per-thread scores into a single vector and returns it
    std::vector<double> getApproxEffectiveResistances() const {
        std::vector<double> result(G.upperNodeIdBound());
        G.parallelForNodes([&](const node u) {
            for (const auto &threadScores : approxEffResistanceGlobal) {
                result[u] += static_cast<double>(threadScores[u]);
            }
            result[u] /= static_cast<double>(numberOfUSTs);
        });

        return result;
    }

    std::vector<double> getDiagonal() { return diagonal; }
    count getRootEccentricity() const { return rootEcc; }
    const Vector &resultVector() const { return result; }
    count computeNumberOfUSTs() const {
        return rootEcc * rootEcc
               * static_cast<count>(
                     std::ceil(std::log(2.0 * static_cast<double>(G.numberOfEdges()) / delta)
                               / (2.0 * epsilon * epsilon)));
    }

    count numberOfUSTs = 0;
    static constexpr count sweeps = 10;
    double getEpsilon() const {
        return epsilon;
    }

    RootStrategy rootStrategy = RootStrategy::MinApproxEcc;

    count solveSingleSystem() {
        ConjugateGradient<CSRMatrix, IdentityPreconditioner> cg(tolerance);
        const auto matrix = CSRMatrix::laplacianMatrix(G);
        cg.setup(matrix);
        Vector rhs(G.upperNodeIdBound());
        rhs[root] = 1.0;
        G.parallelForNodes(
            [&](const node u) { rhs[u] -= 1.0 / static_cast<double>(G.numberOfNodes()); });
        auto status = cg.solve(rhs, result);
        return status.numIters;
    }

    void init();

private:
    // Input parameters
    const Graph &G;
    const double epsilon, delta, tolerance;
    count sampledUSTs = 0;
    node root;
    uint32_t rootEcc;
    Vector result;
    bool didInit = false;

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
    std::vector<std::vector<node>> sequences;

    // Pointers to the parent of the UST, one vector per thread
    std::vector<std::vector<node>> parentGlobal;

    // Index of the parent component of the current component (after the topological order has been
    // determined)
    std::vector<index> biParent;
    // Node within the bionnected component that points to the node in the parent component
    std::vector<node> biAnchor;

    // Topological order of the biconencted components
    std::vector<index> topOrder;

    // Non-normalized approx effective resistance, one per thread.
    std::vector<std::vector<int>> approxEffResistanceGlobal;

    // Timestamps for DFS
    std::vector<std::vector<uint32_t>> tVisitGlobal, tFinishGlobal;

    // Random number generators
    std::vector<pcg32> generators;

    // Parent pointers of the bfs tree
    std::vector<node> bfsParent;

    std::vector<double> diagonal;

    // Creates an entry for the map that contains the bfs-tree edges:
    static constexpr std::pair<std::pair<node, node>, bool> edgeToMapEntry(node u,
                                                                           node v) noexcept {
        return u < v ? std::make_pair(std::make_pair(u, v), true)
                     : std::make_pair(std::make_pair(v, u), false);
    }

    bool isBFSEdge(node u, node v) const { return bfsParent[u] == v || bfsParent[v] == u; }

    // Nodes sequences: Wilson's algorithm runs faster if we start the random walks following a
    // specific sequence of nodes. In this function we compute those sequences.
    void computeNodeSequence();

    void computeBFSTree();
    void sampleUST();
    void dfsUST();
    void aggregateUST();

    void computeDiagonal() {
        solveSingleSystem();
        const auto r = getApproxEffectiveResistances();
        G.parallelForNodes(
            [&](const node u) { diagonal[u] = r[u] - result[root] + 2 * result[u]; });
    }

    node approxMinEccNode();

    // Debugging methods
    void checkBFSTree() const;
    void checkUST() const;
    void checkTwoNodesSequence(const std::vector<node> &sequence) const;
    void checkTimeStamps() const;
};

} // namespace NetworKit

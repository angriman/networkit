/*
 * ForestCentrality.hpp
 *
 *  Created on: 12.02.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#ifndef NETWORKIT_CENTRALITY_FOREST_CENTRALITY
#define NETWORKIT_CENTRALITY_FOREST_CENTRALITY

#include <cmath>
#include <omp.h>
#include <functional>
#include <random>
#include <vector>

#include <ttmath/ttmathbig.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

using high_prec_float = double; //ttmath::Big<//boost::multiprecision::number<boost::multiprecision::cpp_dec_float<2000>>;

// NOTE: Does not support multiple calls of run()
class ForestCentrality final : public Algorithm {

public:
    explicit ForestCentrality(const Graph &G, node root, double epsilon = 0.1);

    void run() override;

    std::vector<count> getNonNormalizedData() const {
        std::vector<count> aggregated(G.upperNodeIdBound());
        G.parallelForNodes([&](const node u) {
            for (const auto &threadScores : approxEffResistanceGlobal)
                aggregated[u] += threadScores[u];
        });
        return aggregated;
    }

    // Aggregates per-thread scores into a single vector and returns it
    std::vector<double> getApproxEffectiveResistance() const {
        std::vector<double> result(G.upperNodeIdBound());
        if (G.isWeighted())
            G.parallelForNodes([&](const node u) {
                high_prec_float tmp = 0;
                for (const auto &threadScores : approxEffResistanceWeightedGlobal)
                    tmp += threadScores[u];
                result[u] = static_cast<double>(tmp / ustsTotalWeight);
            });
        else
            G.parallelForNodes([&](const node u) {
                for (const auto &threadScores : approxEffResistanceGlobal)
                    result[u] += static_cast<double>(threadScores[u]);
                result[u] /= static_cast<double>(numberOfUSTs);
            });

        return result;
    }

    count computeNumberOfUSTs() const {
        return count{4}
               * static_cast<count>(std::ceil(
                   std::log(2.0 * static_cast<double>(G.numberOfEdges()) * G.numberOfNodes())
                   / (2.0 * epsilon * epsilon)));
    }

    count numberOfUSTs = 0;
    double getEpsilon() const { return epsilon; }
    high_prec_float ustsTotalWeight;

private:
    // Input parameters
    const Graph &G;
    node root;
    const double epsilon;

    // Used to mark the status of each node, one vector per thread
    std::vector<std::vector<uint8_t>> statusGlobal;

    // Pointers to the parent of the UST, one vector per thread
    std::vector<std::vector<node>> parentsGlobal;

    // Weights of the edges from each node to their parent in the UST, one vector per thread
    std::vector<std::vector<edgeweight>> weightsToParentGlobal;

    // Nodes sorted by decreasing degree
    std::vector<node> decDegree;

    // Non-normalized approx effective resistance, one per thread.
    std::vector<std::vector<count>> approxEffResistanceGlobal;
    std::vector<std::vector<high_prec_float>> approxEffResistanceWeightedGlobal;

    // Random number generators
    std::vector<std::reference_wrapper<std::mt19937_64>> generators;
    std::vector<std::discrete_distribution<node>> weightedDistr;
    std::vector<std::uniform_int_distribution<count>> uniformDistr;

    void sampleUSTsUnweighted();
    void sampleUSTsWeighted();
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_FOREST_CENTRALITY

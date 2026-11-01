// networkit-format
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/centrality/HITS.hpp>

namespace NetworKit {

HITS::HITS(const Graph &G, double tol, count maxIterations)
    : Centrality(G), tol(tol), maxIterations(maxIterations) {
    if (tol <= 0) {
        throw std::invalid_argument("Error: tol must be positive.");
    }
    if (maxIterations == 0) {
        throw std::invalid_argument("Error: maxIterations must be positive.");
    }
}

void HITS::run() {
    Aux::SignalHandler handler;
    const auto n = G.numberOfNodes();
    const auto z = G.upperNodeIdBound();

    scoreData.assign(z, 0.);
    hubData.assign(z, 0.);
    iterations = 0;

    if (n == 0) {
        hasRun = true;
        return;
    }

    const double initialScore = 1. / std::sqrt(static_cast<double>(n));
    G.parallelForNodes([&](node u) -> void {
        scoreData[u] = initialScore;
        hubData[u] = initialScore;
    });
    std::vector<double> nextAuthority(z), nextHub(z);

    auto normalize = [&](std::vector<double> &values) -> void {
        const double norm = std::sqrt(
            G.parallelSumForNodes([&](node u) -> double { return values[u] * values[u]; }));
        if (norm > 0) {
            G.parallelForNodes([&](node u) -> void { values[u] /= norm; });
        }
    };

    double distance;
    do {
        handler.assureRunning();

        G.balancedParallelForNodes([&](node u) -> void {
            nextAuthority[u] = 0.;
            G.forInEdgesOf(
                u, [&](node, node v, edgeweight w) -> void { nextAuthority[u] += w * hubData[v]; });
        });
        normalize(nextAuthority);

        G.balancedParallelForNodes([&](node u) -> void {
            nextHub[u] = 0.;
            G.forEdgesOf(
                u, [&](node, node v, edgeweight w) -> void { nextHub[u] += w * nextAuthority[v]; });
        });
        normalize(nextHub);

        distance = std::sqrt(G.parallelSumForNodes([&](node u) -> double {
            const double authorityDiff = nextAuthority[u] - scoreData[u];
            const double hubDiff = nextHub[u] - hubData[u];
            return authorityDiff * authorityDiff + hubDiff * hubDiff;
        }));

        std::swap(scoreData, nextAuthority);
        std::swap(hubData, nextHub);
        ++iterations;
    } while (distance > tol && iterations < maxIterations);

    hasRun = true;
}

} /* namespace NetworKit */

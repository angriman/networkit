// networkit-format
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/community/StochasticBlockmodelInference.hpp>

namespace NetworKit {

namespace {

// Moves that gain less are rounding noise, and accepting them could cycle.
constexpr double minGain = 1e-9;

double xLogX(double x) {
    return x > 0. ? x * std::log(x) : 0.;
}

/*
 * Block assignment with the statistics of the likelihood. Writing m(r, s) for the edge weight
 * between blocks r and s (twice the internal weight if r = s), kappa(r) for the total degree of r,
 * and size(r) for its number of nodes, the likelihood of Karrer and Newman simplifies to
 *
 *     L = sum_{r,s} m(r, s) log m(r, s) - 2 sum_r kappa(r) log kappa(r)
 *
 * and, without degree correction, to the same expression with log size(r) in the second term.
 * Moving a node therefore only changes the entries of m next to its source and target blocks.
 */
class BlockState {
public:
    BlockState(const Graph &G, count k, bool degreeCorrected)
        : G(G), k(k), degreeCorrected(degreeCorrected), block(G.upperNodeIdBound(), none),
          weights(k * k), degrees(k), sizes(k), neighborWeight(k, 0.) {
        adjacentBlocks.reserve(k);
    }

    void assign(const std::vector<index> &newBlock) {
        block = newBlock;
        std::fill(weights.begin(), weights.end(), 0.);
        std::fill(degrees.begin(), degrees.end(), 0.);
        std::fill(sizes.begin(), sizes.end(), count{0});

        G.forNodes([&](node u) -> void {
            ++sizes[block[u]];
            degrees[block[u]] += G.weightedDegree(u, true);
            G.forEdgesOf(u, [&](node, node v, edgeweight weight, edgeid) -> void {
                // Every edge is visited from both ends, and a self-loop counts twice in m.
                weightBetween(block[u], block[v]) += u == v ? 2. * weight : weight;
            });
        });
    }

    double logLikelihood() const {
        double result = 0.;
        for (const auto weight : weights) {
            result += xLogX(weight);
        }
        for (index r = 0; r < k; ++r) {
            result -= 2. * penalty(degrees[r], sizes[r]);
        }
        return result;
    }

    index blockOf(node u) const { return block[u]; }

    const std::vector<index> &blocks() const { return block; }

    // Must be called before moveGain() and move() for the same node.
    void selectNode(node u) {
        for (const auto r : adjacentBlocks) {
            neighborWeight[r] = 0.;
        }
        adjacentBlocks.clear();
        selected = u;
        selfLoopWeight = 0.;
        G.forEdgesOf(u, [&](node, node v, edgeweight weight, edgeid) -> void {
            if (weight == 0.) {
                return;
            }
            if (u == v) {
                selfLoopWeight = 2. * weight;
                return;
            }
            const index r = block[v];
            if (neighborWeight[r] == 0.) {
                adjacentBlocks.push_back(r);
            }
            neighborWeight[r] += weight;
        });
    }

    // Returns the best block other than the current one of the selected node, and the gain.
    std::pair<index, double> bestMove(count numberOfBlocks) const {
        index bestTarget = none;
        double bestGain = -std::numeric_limits<double>::infinity();
        for (index target = 0; target < numberOfBlocks; ++target) {
            if (target == block[selected]) {
                continue;
            }
            const double gain = moveGain(target);
            if (gain > bestGain) {
                bestGain = gain;
                bestTarget = target;
            }
        }
        return {bestTarget, bestGain};
    }

    double moveGain(index target) const {
        const index source = block[selected];
        const double toSource = neighborWeight[source];
        const double toTarget = neighborWeight[target];

        double gain = 0.;
        for (const auto t : adjacentBlocks) {
            if (t != source && t != target) {
                // Each off-diagonal entry appears twice in the sum over ordered pairs.
                gain += 2.
                        * (change(weightBetween(source, t), -neighborWeight[t])
                           + change(weightBetween(target, t), neighborWeight[t]));
            }
        }
        gain += 2. * change(weightBetween(source, target), toSource - toTarget);
        gain += change(weightBetween(source, source), -2. * toSource - selfLoopWeight);
        gain += change(weightBetween(target, target), 2. * toTarget + selfLoopWeight);

        const double degree = G.weightedDegree(selected, true);
        gain -=
            2.
            * (penalty(degrees[source] - degree, sizes[source] - 1)
               + penalty(degrees[target] + degree, sizes[target] + 1)
               - penalty(degrees[source], sizes[source]) - penalty(degrees[target], sizes[target]));
        return gain;
    }

    void move(index target) {
        const index source = block[selected];
        const double toSource = neighborWeight[source];
        const double toTarget = neighborWeight[target];

        for (const auto t : adjacentBlocks) {
            if (t != source && t != target) {
                addWeight(source, t, -neighborWeight[t]);
                addWeight(target, t, neighborWeight[t]);
            }
        }
        addWeight(source, target, toSource - toTarget);
        weightBetween(source, source) -= 2. * toSource + selfLoopWeight;
        weightBetween(target, target) += 2. * toTarget + selfLoopWeight;

        const double degree = G.weightedDegree(selected, true);
        degrees[source] -= degree;
        degrees[target] += degree;
        --sizes[source];
        ++sizes[target];
        block[selected] = target;
    }

private:
    double &weightBetween(index r, index s) { return weights[r * k + s]; }

    double weightBetween(index r, index s) const { return weights[r * k + s]; }

    void addWeight(index r, index s, double delta) {
        weightBetween(r, s) += delta;
        weightBetween(s, r) += delta;
    }

    static double change(double value, double delta) { return xLogX(value + delta) - xLogX(value); }

    double penalty(double degree, count size) const {
        if (degreeCorrected) {
            return xLogX(degree);
        }
        return degree > 0. ? degree * std::log(static_cast<double>(size)) : 0.;
    }

    const Graph &G;
    count k;
    bool degreeCorrected;

    std::vector<index> block;
    std::vector<double> weights;
    std::vector<double> degrees;
    std::vector<count> sizes;

    node selected = none;
    double selfLoopWeight = 0.;
    std::vector<double> neighborWeight;
    std::vector<index> adjacentBlocks;
};

// Moves every node to its best block if that improves the likelihood.
bool greedyPass(BlockState &state, const std::vector<node> &order, count k) {
    bool improved = false;
    for (const auto u : order) {
        state.selectNode(u);
        const auto [target, gain] = state.bestMove(k);
        if (gain > minGain) {
            state.move(target);
            improved = true;
        }
    }
    return improved;
}

} // namespace

StochasticBlockmodelInference::StochasticBlockmodelInference(const Graph &G, count k,
                                                             bool degreeCorrected, count restarts,
                                                             count maxSweeps)
    : CommunityDetectionAlgorithm(G), k(k), degreeCorrected(degreeCorrected), restarts(restarts),
      maxSweeps(maxSweeps) {
    if (k == 0 || k > G.numberOfNodes()) {
        throw std::invalid_argument(
            "StochasticBlockmodelInference: k must be between 1 and the number of nodes.");
    }
    if (restarts == 0 || maxSweeps == 0) {
        throw std::invalid_argument(
            "StochasticBlockmodelInference: restarts and maxSweeps must be at least 1.");
    }
    G.forEdges([](node, node, edgeweight weight) -> void {
        if (weight < 0.) {
            throw std::invalid_argument(
                "StochasticBlockmodelInference: edge weights must not be negative.");
        }
    });
}

void StochasticBlockmodelInference::run() {
    std::vector<node> nodes;
    nodes.reserve(G->numberOfNodes());
    G->forNodes([&](node u) -> void { nodes.push_back(u); });

    double bestLikelihood = -std::numeric_limits<double>::infinity();
    std::vector<index> bestBlocks;

#pragma omp parallel
    {
        BlockState state(*G, k, degreeCorrected);
        std::vector<node> order(nodes);
        std::vector<index> start(G->upperNodeIdBound(), none);

#pragma omp for schedule(dynamic)
        for (omp_index restart = 0; restart < static_cast<omp_index>(restarts); ++restart) {
            for (const auto u : nodes) {
                start[u] = Aux::Random::integer(k - 1);
            }
            state.assign(start);

            for (count sweep = 0; sweep < maxSweeps; ++sweep) {
                std::shuffle(order.begin(), order.end(), Aux::Random::getURNG());
                if (!greedyPass(state, order, k)) {
                    break;
                }
            }

            // Recomputing discards the rounding errors of the incremental updates.
            state.assign(state.blocks());
            const double likelihood = state.logLikelihood();
#pragma omp critical
            {
                if (likelihood > bestLikelihood) {
                    bestLikelihood = likelihood;
                    bestBlocks = state.blocks();
                }
            }
        }
    }

    result = Partition(G->upperNodeIdBound());
    result.setUpperBound(k);
    G->forNodes([&](node u) -> void { result[u] = bestBlocks[u]; });
    logLikelihood = bestLikelihood;
    hasRun = true;
}

double StochasticBlockmodelInference::getLogLikelihood() const {
    assureFinished();
    return logLikelihood;
}

} // namespace NetworKit

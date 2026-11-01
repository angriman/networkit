// networkit-format
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/community/StochasticBlockmodelInference.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

namespace {

// Evaluates the likelihood of Karrer and Newman from its original form.
double likelihoodFromDefinition(const Graph &G, const std::vector<index> &block, count k,
                                bool degreeCorrected) {
    std::vector<std::vector<double>> m(k, std::vector<double>(k, 0.));
    std::vector<double> kappa(k, 0.);
    std::vector<double> sizes(k, 0.);
    G.forNodes([&](node u) -> void { ++sizes[block[u]]; });
    G.forEdges([&](node u, node v, edgeweight weight) -> void {
        m[block[u]][block[v]] += weight;
        m[block[v]][block[u]] += weight;
        kappa[block[u]] += weight;
        kappa[block[v]] += weight;
    });

    double result = 0.;
    for (index r = 0; r < k; ++r) {
        for (index s = 0; s < k; ++s) {
            if (m[r][s] > 0.) {
                const double expected = degreeCorrected ? kappa[r] * kappa[s] : sizes[r] * sizes[s];
                result += m[r][s] * std::log(m[r][s] / expected);
            }
        }
    }
    return result;
}

double bestLikelihood(const Graph &G, count k, bool degreeCorrected) {
    std::vector<node> nodes;
    G.forNodes([&](node u) -> void { nodes.push_back(u); });
    std::vector<index> block(G.upperNodeIdBound(), 0);
    double best = -std::numeric_limits<double>::infinity();
    while (true) {
        best = std::max(best, likelihoodFromDefinition(G, block, k, degreeCorrected));
        index position = 0;
        while (position < nodes.size() && ++block[nodes[position]] == k) {
            block[nodes[position]] = 0;
            ++position;
        }
        if (position == nodes.size()) {
            return best;
        }
    }
}

} // namespace

TEST(StochasticBlockmodelInferenceGTest, testInvalidArguments) {
    Graph G(4);
    G.addEdge(0, 1);
    EXPECT_THROW(StochasticBlockmodelInference(G, 0), std::invalid_argument);
    EXPECT_THROW(StochasticBlockmodelInference(G, 5), std::invalid_argument);
    EXPECT_THROW(StochasticBlockmodelInference(G, 2, true, 0), std::invalid_argument);
    EXPECT_THROW(StochasticBlockmodelInference(G, 2, true, 1, 0), std::invalid_argument);

    Graph directed(2, false, true);
    EXPECT_THROW(StochasticBlockmodelInference(directed, 1), std::runtime_error);

    Graph negative(2, true);
    negative.addEdge(0, 1, -1.);
    EXPECT_THROW(StochasticBlockmodelInference(negative, 1), std::invalid_argument);
}

TEST(StochasticBlockmodelInferenceGTest, testNotRun) {
    Graph G(2);
    StochasticBlockmodelInference algorithm(G, 1);
    EXPECT_THROW(algorithm.getLogLikelihood(), std::runtime_error);
    EXPECT_THROW(algorithm.getPartition(), std::runtime_error);
}

TEST(StochasticBlockmodelInferenceGTest, testAgainstExhaustiveSearch) {
    Aux::Random::setSeed(42, true);
    for (const bool weighted : {false, true}) {
        for (const bool degreeCorrected : {true, false}) {
            for (count trial = 0; trial < 5; ++trial) {
                Graph G = ErdosRenyiGenerator(8, 0.4).generate();
                G = Graph(G, weighted, false);
                if (weighted) {
                    G.forEdges([&](node u, node v) -> void {
                        G.setWeight(u, v, static_cast<edgeweight>(Aux::Random::integer(1, 5)));
                    });
                }
                for (const count k : {2, 3}) {
                    StochasticBlockmodelInference algorithm(G, k, degreeCorrected, 100);
                    algorithm.run();
                    EXPECT_NEAR(algorithm.getLogLikelihood(), bestLikelihood(G, k, degreeCorrected),
                                1e-9);
                }
            }
        }
    }
}

TEST(StochasticBlockmodelInferenceGTest, testKarateClub) {
    Aux::Random::setSeed(42, true);
    const Graph G = METISGraphReader{}.read("input/karate.graph");

    const auto blockOfFirstNode = [&](bool degreeCorrected) -> std::vector<node> {
        StochasticBlockmodelInference algorithm(G, 2, degreeCorrected, 200);
        algorithm.run();
        const auto &partition = algorithm.getPartition();
        std::vector<node> result;
        G.forNodes([&](node u) -> void {
            if (partition[u] == partition[0]) {
                result.push_back(u);
            }
        });
        return result;
    };

    // As reported by Karrer and Newman: the plain model separates the hubs from the rest, while
    // the degree-corrected one finds the two factions, up to nodes 8 and 9.
    EXPECT_EQ(blockOfFirstNode(false), (std::vector<node>{0, 1, 2, 32, 33}));
    EXPECT_EQ(blockOfFirstNode(true),
              (std::vector<node>{0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 16, 17, 19, 21}));
}

TEST(StochasticBlockmodelInferenceGTest, testBipartiteGraph) {
    Aux::Random::setSeed(42, true);
    // There are no communities here in the usual sense, as all edges run between the two sides.
    const count n = 40;
    Graph G(n);
    for (node u = 0; u < n / 2; ++u) {
        for (node v = n / 2; v < n; ++v) {
            if (Aux::Random::probability() < 0.3) {
                G.addEdge(u, v);
            }
        }
    }

    StochasticBlockmodelInference algorithm(G, 2);
    algorithm.run();
    const auto &partition = algorithm.getPartition();
    G.forNodes([&](node u) -> void { EXPECT_EQ(partition[u] == partition[0], u < n / 2); });
}

TEST(StochasticBlockmodelInferenceGTest, testLikelihoodOfGeneralGraph) {
    Aux::Random::setSeed(42, true);
    Graph G(12, true);
    G.addEdge(0, 1, 2.5);
    G.addEdge(1, 2, 1.5);
    G.addEdge(0, 2, 3.);
    G.addEdge(2, 2, 2.);
    G.addEdge(3, 4, 0.);
    G.addEdge(4, 5, 4.);
    G.addEdge(5, 6, 1.);
    G.addEdge(6, 6, 0.5);
    G.addEdge(6, 7, 2.);
    G.addEdge(8, 9, 3.);
    G.addEdge(2, 3, 0.25);
    G.removeNode(10);

    for (const bool degreeCorrected : {true, false}) {
        for (const count k : {1, 2, 3, 4}) {
            StochasticBlockmodelInference algorithm(G, k, degreeCorrected);
            algorithm.run();
            const auto &partition = algorithm.getPartition();
            EXPECT_EQ(partition[10], none);
            G.forNodes([&](node u) -> void { EXPECT_LT(partition[u], k); });
            EXPECT_NEAR(algorithm.getLogLikelihood(),
                        likelihoodFromDefinition(G, partition.getVector(), k, degreeCorrected),
                        1e-9);
        }
    }
}

} // namespace NetworKit

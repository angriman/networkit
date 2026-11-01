// networkit-format
/*
 * DegreeSequenceGeneratorGTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <algorithm>
#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ConfigurationModel.hpp>
#include <networkit/generators/EdgeSwitchingMarkovChainGenerator.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/generators/PowerlawDegreeSequence.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class DegreeSequenceGeneratorGTest : public testing::Test {};

TEST_F(DegreeSequenceGeneratorGTest,
       testConfigurationModelGeneratorWithEdgeSwitchingOnRealSequence) {
    METISGraphReader reader;
    auto graphs = {"input/jazz.graph", "input/lesmis.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.upperNodeIdBound();
        std::vector<count> sequence(n);
        G.forNodes([&](node u) { sequence[u] = G.degree(u); });

        bool skipTest = false;
        EdgeSwitchingMarkovChainGenerator gen(sequence, skipTest);
        Graph G2 = gen.generate();

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        std::vector<count> testSequence(n);
        G2.forNodes([&](node u) { testSequence[u] = G2.degree(u); });
        Aux::Parallel::sort(testSequence.begin(), testSequence.end(), std::greater<count>());
        Aux::Parallel::sort(sequence.begin(), sequence.end(), std::greater<count>());

        for (index i = 0; i < n; ++i) {
            EXPECT_EQ(sequence[i], testSequence[i]);
        }
    }
}

TEST_F(DegreeSequenceGeneratorGTest,
       testConfigurationModelGeneratorWithRejectionSamplingOnRealSequence) {
    METISGraphReader reader;
    Aux::Random::setSeed(42, false);
    auto graphs = {"input/jazz.graph", "input/lesmis.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.numberOfNodes();
        std::vector<count> sequence(n);
        G.forNodes([&](node u) { sequence[u] = G.degree(u); });

        ConfigurationModel gen(sequence);
        Graph G2 = gen.generate();

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        std::vector<count> testSequence(n);
        G2.forNodes([&](node u) { testSequence[u] = G2.degree(u); });
        Aux::Parallel::sort(testSequence.begin(), testSequence.end(), std::greater<count>());
        Aux::Parallel::sort(sequence.begin(), sequence.end(), std::greater<count>());

        for (index i = 0; i < n; ++i) {
            EXPECT_EQ(sequence[i], testSequence[i]);
        }
    }
}

TEST_F(DegreeSequenceGeneratorGTest, testHavelHakimiGeneratorOnRandomSequence) {
    count n = 400;
    std::vector<count> sequence(n);
    bool realizable = false;

    std::mt19937_64 prng(1);
    std::uniform_int_distribution<count> deg_distr(0, n / 10 - 1);
    do {
        // fill sequence with random values (this is not power-law, of course!)
        std::generate(sequence.begin(), sequence.end(), [&] { return deg_distr(prng); });

        // check if sequence is realizable
        HavelHakimiGenerator hhgen(sequence);
        realizable = hhgen.isRealizable();

        if (realizable) {
            Graph G = hhgen.generate();
            EXPECT_TRUE(G.checkConsistency());
            count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
            EXPECT_EQ(volume, 2 * G.numberOfEdges());
        }
    } while (!realizable);
}

TEST_F(DegreeSequenceGeneratorGTest, testHavelHakimiGeneratorOnRealSequence) {
    METISGraphReader reader;
    auto graphs = {"input/jazz.graph", "input/lesmis.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.numberOfNodes();
        std::vector<count> sequence(n);
        G.forNodes([&](node u) { sequence[u] = G.degree(u); });

        HavelHakimiGenerator hhgen(sequence);
        Graph G2 = hhgen.generate();
        EXPECT_TRUE(G.checkConsistency());

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        if (volume < 50000) {
            std::vector<count> testSequence(n);
            G2.forNodes([&](node u) { testSequence[u] = G2.degree(u); });

            for (index i = 0; i < n; ++i) {
                EXPECT_EQ(sequence[i], testSequence[i]);
            }
        }
    }
}

TEST_F(DegreeSequenceGeneratorGTest, testHavelHakimiGeneratorOnUnrealizableSequence) {
    std::vector<count> seq = {20, 10, 2, 2, 2, 2, 2, 2, 2, 2, 2};

    HavelHakimiGenerator hhgen(seq);
    EXPECT_THROW(hhgen.generate(), std::runtime_error);

    hhgen = HavelHakimiGenerator(seq, true);
    Graph G = hhgen.generate();

    G.forNodes([&](node u) { EXPECT_EQ(std::min<count>(seq[u], 10), G.degree(u)); });
}

TEST_F(DegreeSequenceGeneratorGTest, testStaticDegreeSequenceGenerator) {
    auto test_known = [](const std::vector<count> &seq, bool result) {
        HavelHakimiGenerator gen(seq, true);
        ASSERT_EQ(gen.isRealizable(), result);
    };

    // some contrived sequence
    test_known({1}, false);
    test_known({1, 1}, true);
    test_known({2, 2, 2}, true);
    test_known({1, 2, 1}, true);
    test_known({1, 3, 1}, false);
    test_known({1, 1, 3, 1}, true);

    // some random sequences
    std::mt19937_64 prng(1);
    std::uniform_int_distribution<node> distr_num_nodes(1, 50);

    unsigned numRealized = 0;
    for (int iter = 0; iter < 100; ++iter) {
        const auto n = distr_num_nodes(prng);
        std::vector<count> seq(n);
        std::generate(seq.begin(), seq.end(),
                      [&] { return std::uniform_int_distribution<count>{0, n - 1}(prng); });

        HavelHakimiGenerator gen(seq, true);
        const auto isRealizable = gen.isRealizable();

        const auto G = gen.generate();
        const auto nodes = G.nodeRange();
        ASSERT_EQ(G.numberOfNodes(), n);

        const auto didRealize =
            std::all_of(nodes.begin(), nodes.end(), [&](node u) { return G.degree(u) == seq[u]; });

        ASSERT_EQ(isRealizable, didRealize);
        numRealized += didRealize;
    }

    ASSERT_GT(numRealized, 10);
    ASSERT_LT(numRealized, 90);
}

TEST_F(DegreeSequenceGeneratorGTest, testConfigurationModelGeneratorOnRealSequence) {
    METISGraphReader reader;
    auto graphs = {"input/jazz.graph", "input/lesmis.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.numberOfNodes();
        std::vector<count> sequence(n);
        G.forNodes([&](node u) { sequence[u] = G.degree(u); });

        bool skipTest = false;
        EdgeSwitchingMarkovChainGenerator gen(sequence, skipTest);
        Graph G2 = gen.generate();

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        if (volume < 50000) {
            std::vector<count> testSequence(n);
            G2.forNodes([&](node u) { testSequence[u] = G2.degree(u); });
            Aux::Parallel::sort(testSequence.begin(), testSequence.end(), std::greater<count>());
            Aux::Parallel::sort(sequence.begin(), sequence.end(), std::greater<count>());

            for (index i = 0; i < n; ++i) {
                EXPECT_EQ(sequence[i], testSequence[i]);
            }
        }
    }
}

TEST_F(DegreeSequenceGeneratorGTest, testPowerLawDegreeSequenceFromDegreeSequence) {
    std::vector<double> inputVector = {2, 4, 3, 2, 1, 6};
    double sum = 0.0;
    for (auto it = inputVector.begin(); it != inputVector.end(); ++it) {
        sum += *it;
    }
    auto avg = sum / inputVector.size();

    PowerlawDegreeSequence PLDS(inputVector);
    PLDS.setGammaFromAverageDegree(avg);
    PLDS.run();
    EXPECT_NEAR(PLDS.getExpectedAverageDegree(), avg, 1.0);
    EXPECT_EQ(PLDS.getMinimumDegree(), 1);
    EXPECT_EQ(PLDS.getMaximumDegree(), 6);
}

TEST_F(DegreeSequenceGeneratorGTest, testPowerLawDegreeSequenceFromGraph) {
    METISGraphReader reader;
    auto graphpath = "input/jazz.graph";
    Graph G = reader.read(graphpath);

    PowerlawDegreeSequence PLDS(G);
    PLDS.run();
    EXPECT_NEAR(PLDS.getExpectedAverageDegree(), 19, 1);
    EXPECT_EQ(PLDS.getMinimumDegree(), 1);
    EXPECT_EQ(PLDS.getMaximumDegree(), 100);
}

} /* namespace NetworKit */

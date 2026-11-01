// networkit-format
/*
 * LFRGeneratorGTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <fstream>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/LFRGenerator.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class LFRGeneratorGTest : public testing::Test {
protected:
    template <class Type>
    std::vector<Type> readVector(std::string_view path) const {
        std::ifstream inputFile(path.data());
        Type cur;
        std::vector<Type> data;

        while (inputFile >> cur)
            data.push_back(cur);

        return data;
    }
};

TEST_F(LFRGeneratorGTest, testLFRGenerator) {
    Aux::Random::setSeed(42, false);
    count n = 500;
    LFRGenerator gen(n);
    gen.generatePowerlawDegreeSequence(20, 50, -2);
    gen.generatePowerlawCommunitySizeSequence(10, 50, -1);
    gen.setMu(0.5);
    gen.run();
    Graph G1 = gen.getMoveGraph();
    gen.run(); // should rewire the edges but nothing else
    Graph G2 = gen.getMoveGraph();
    EXPECT_EQ(n, G1.numberOfNodes());
    EXPECT_EQ(n, G2.numberOfNodes());
    EXPECT_EQ(G1.numberOfEdges(), G2.numberOfEdges());
}

TEST_F(LFRGeneratorGTest, testLFRGeneratorImpossibleSequence) {
    LFRGenerator gen(100);
    gen.generatePowerlawDegreeSequence(10, 11, -2);
    EXPECT_ANY_THROW(gen.generatePowerlawCommunitySizeSequence(9, 8, -3));
    gen.setMu(0.5);
    EXPECT_THROW(gen.run(), std::runtime_error);
    EXPECT_THROW(gen.getMoveGraph(), std::runtime_error);
}

TEST_F(LFRGeneratorGTest, testLFRGeneratorWithRealData) {
    const auto degreeSequence = readVector<count>("input/testLFR/testLFRDegSeq.txt"),
               partition = readVector<count>("input/testLFR/testLFRPartition.txt");

    const auto mu = readVector<double>("input/testLFR/testLFRMu.txt");
    Partition C(partition.size());
    C.setUpperBound(20);
    for (node u = 0; u < partition.size(); ++u) {
        C[u] = partition[u];
    }
    LFRGenerator gen(degreeSequence.size());
    gen.setDegreeSequence(degreeSequence);
    gen.setPartition(C);
    gen.setMu(mu);
    gen.run();
    Graph G = gen.getGraph();
    G.parallelForNodes([&](node u) { EXPECT_EQ(G.degree(u), degreeSequence[u]); });
    EXPECT_EQ(C.numberOfSubsets(), gen.getPartition().numberOfSubsets());
}

TEST_F(LFRGeneratorGTest, testLFRGeneratorWithBinomialDistribution) {
    Aux::Random::setSeed(42, false);
    count n = 500;
    LFRGenerator gen(n);
    gen.generatePowerlawDegreeSequence(20, 50, -2);
    gen.generatePowerlawCommunitySizeSequence(10, 50, -1);
    gen.setMuWithBinomialDistribution(0.5);
    gen.run();
    Graph G1 = gen.getMoveGraph();
    gen.run(); // should rewire the edges but nothing else
    Graph G2 = gen.getMoveGraph();
    EXPECT_EQ(n, G1.numberOfNodes());
    EXPECT_EQ(n, G2.numberOfNodes());
    EXPECT_EQ(G1.numberOfEdges(), G2.numberOfEdges());
}

} /* namespace NetworKit */

// networkit-format
/*
 * ChungLuGeneratorGTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ChungLuGenerator.hpp>
#include <networkit/generators/ChungLuGeneratorAlamEtAl.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class ChungLuGeneratorGTest : public testing::Test {};

TEST_F(ChungLuGeneratorGTest, testChungLuGenerator) {
    count n = 400;
    count maxDegree = n / 8;
    std::vector<count> sequence(n);
    count expVolume = 0;
    count actualVolume = 0;

    // fill sequence with random values (this is not power-law, of course!)
    for (index i = 0; i < n; ++i) {
        sequence[i] = rand() % maxDegree;
        expVolume += sequence[i];
    }

    ChungLuGenerator gen(sequence);
    Graph G = gen.generate();
    EXPECT_TRUE(G.checkConsistency());

    EXPECT_EQ(n, G.numberOfNodes());
    G.forNodes([&](node v) { actualVolume += G.degree(v); });

    INFO("expected volume: ", expVolume, ", actual volume: ", actualVolume);
}

TEST_F(ChungLuGeneratorGTest, testChungLuGeneratorDegreeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++) {
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
    }
    ChungLuGenerator generator(vec);
    Graph G = generator.generate();
    /* We check to see if the actual degrees of our nodes vary too much from the expected ones.
     * However, we need to sort the expected degrees first, since the algorithm does this as well
     * and the nodes with the highest degrees are added first. */
    Aux::Parallel::sort(vec.begin(), vec.end(), [](count a, count b) { return a > b; });
    /* Check if node degree is more than 50% off from the expected degree of that node. */
    // TODO Should we be looking for something better than a 50% range here?
    G.parallelForNodes([&](node v) { EXPECT_NEAR(G.degree(v), vec[v], (0.5 * maxDegree)); });
}

TEST_F(ChungLuGeneratorGTest, testChungLuGeneratorVolumeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    count expectedVolume = 0;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++) {
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
        expectedVolume += grad;
    }
    ChungLuGenerator generator(vec);
    Graph G = generator.generate();
    /* Check if volume is more than 10% off from the expected volume. */
    // TODO Is a 20% offset here sufficient? */
    EXPECT_NEAR(G.numberOfEdges() * 2, expectedVolume, 0.2 * expectedVolume);
}

TEST_F(ChungLuGeneratorGTest, testChungLuGeneratorAlamEtAl) {
    count n = 400;
    count maxDegree = n / 8;
    std::vector<count> sequence(n);
    count expVolume = 0;
    count actualVolume = 0;

    // fill sequence with random values (this is not power-law, of course!)
    for (index i = 0; i < n; ++i) {
        sequence[i] = rand() % maxDegree;
        expVolume += sequence[i];
    }

    ChungLuGeneratorAlamEtAl gen(sequence);
    Graph G = gen.generate();
    EXPECT_TRUE(G.checkConsistency());

    EXPECT_EQ(n, G.numberOfNodes());
    G.forNodes([&](node v) { actualVolume += G.degree(v); });

    INFO("expected volume: ", expVolume, ", actual volume: ", actualVolume);
}

TEST_F(ChungLuGeneratorGTest, testChungLuGeneratorAlamEtAlDegreeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++) {
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
    }
    ChungLuGeneratorAlamEtAl generator(vec);
    Graph G = generator.generate();
    /* We check to see if the actual degrees of our nodes vary too much from the expected ones.
     * However, we need to sort the expected degrees first, since the algorithm does this as well
     * and the nodes with the highest degrees are added first. */
    Aux::Parallel::sort(vec.begin(), vec.end(), [](count a, count b) { return a < b; });
    /* Check if node degree is more than 50% off from the expected degree of that node. */
    // TODO Should we be looking for something better than a 50% range here?
    G.parallelForNodes([&](node v) { EXPECT_NEAR(G.degree(v), vec[v], (0.5 * maxDegree)); });
}

TEST_F(ChungLuGeneratorGTest, testChungLuGeneratorAlamEtAlVolumeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    count expectedVolume = 0;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++) {
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
        expectedVolume += grad;
    }
    ChungLuGeneratorAlamEtAl generator(vec);
    Graph G = generator.generate();
    /* Check if volume is more than 10% off from the expected volume. */
    // TODO Is a 20% offset here sufficient? */
    EXPECT_NEAR(G.numberOfEdges() * 2, expectedVolume, 0.2 * expectedVolume);
}

} /* namespace NetworKit */

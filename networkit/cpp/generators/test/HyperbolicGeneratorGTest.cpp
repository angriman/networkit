// networkit-format
/*
 * HyperbolicGeneratorGTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/dynamics/GraphUpdater.hpp>
#include <networkit/generators/DynamicHyperbolicGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/viz/PostscriptWriter.hpp>

namespace NetworKit {

class HyperbolicGeneratorGTest : public testing::Test {
public:
    std::vector<double> getAngles(const DynamicHyperbolicGenerator &dynGen) {
        return dynGen.angles;
    }

    std::vector<double> getRadii(const DynamicHyperbolicGenerator &dynGen) { return dynGen.radii; }
};

/**
 * Testing the dynamic hyperbolic generator with fixed parameters and changing node positions
 */
TEST_F(HyperbolicGeneratorGTest, testDynamicHyperbolicGeneratorOnMovedNodes) {
    Aux::Random::setSeed(0, false);

    // set up dynamic parameters
    int nSteps = 20;
    const count n = 500;
    const double k = 6;
    const double alpha = 1;
    // const double exp = 2*alpha+1;
    const double T = 0;
    const double R = HyperbolicSpace::getTargetRadius(n, n * k / 2, alpha, T);

    double movedShare = 1;
    double moveDistance = 0.1;

    // set up initial node positions
    vector<double> angles(n, -1);
    vector<double> radii(n, -1);
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);
    DynamicHyperbolicGenerator dynGen(angles, radii, R, alpha, T, movedShare, moveDistance);

    // generate starting graph
    Graph G = HyperbolicGenerator().generate(angles, radii, R);
    count initialEdgeCount = G.numberOfEdges();
    count expected = n * HyperbolicSpace::getExpectedDegree(n, alpha, R) * 0.5;
    EXPECT_NEAR(initialEdgeCount, expected, expected / 5);
    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;

    for (int i = 0; i < nSteps; i++) {
        // move nodes and generate stream of affected edges
        stream = dynGen.generate(1);
        DEBUG("Edges: ", G.numberOfEdges());
        for (auto event : stream) {
            EXPECT_TRUE(event.type == GraphEvent::EDGE_REMOVAL
                        || event.type == GraphEvent::EDGE_ADDITION
                        || event.type == GraphEvent::TIME_STEP);
            if (event.type == GraphEvent::EDGE_REMOVAL) {
                EXPECT_TRUE(G.hasEdge(event.u, event.v));
            }
            // only present nodes can be affected, no new nodes are introduced
            if (event.type != GraphEvent::TIME_STEP) {
                EXPECT_LT(event.u, G.upperNodeIdBound());
            }
        }
        gu.update(stream);
        EXPECT_TRUE(G.checkConsistency());
    }

    // update moved nodes
    angles = getAngles(dynGen);
    radii = getRadii(dynGen);
    Graph comparison = HyperbolicGenerator().generate(angles, radii, R);
    EXPECT_EQ(G.numberOfEdges(), comparison.numberOfEdges());

    // heuristic criterion: Number of edges may change, but should not change much
    EXPECT_NEAR(G.numberOfEdges(), initialEdgeCount, initialEdgeCount / 5);
}

/**
 * creates a series of pictures visualizing the effect of the dynamic hyperbolic generator
 */
TEST_F(HyperbolicGeneratorGTest, testDynamicHyperbolicVisualization) {
    count n = 300;
    count nSteps = 20;

    const double k = 6;
    const double alpha = 1;
    // const double exp = 2*alpha+1;
    const double T = 0;
    const double R = HyperbolicSpace::getTargetRadius(n, n * k / 2, alpha, T);

    double movedShare = 0.2;
    double moveDistance = 1;
    vector<double> angles(n);
    vector<double> radii(n);

    HyperbolicSpace::fillPoints(angles, radii, R, alpha);

    DynamicHyperbolicGenerator dynGen(angles, radii, R, alpha, T, movedShare, moveDistance);
    Graph G = dynGen.getGraph();

    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;

    for (index i = 0; i < nSteps; i++) {
        stream = dynGen.generate(1);
        DEBUG("Edges: ", G.numberOfEdges());
        for (auto event : stream) {
            EXPECT_TRUE(event.type == GraphEvent::EDGE_REMOVAL
                        || event.type == GraphEvent::EDGE_ADDITION
                        || event.type == GraphEvent::TIME_STEP);
        }
        gu.update(stream);

        auto coordinates = dynGen.getCoordinates();

#ifndef NETWORKIT_RELEASE_LOGGING
        {
            std::stringstream ss;
            ss << "output/hyperbolic-" << std::setw(4) << std::setfill('0') << i << ".eps";

            PostscriptWriter psWriter(true);
            psWriter.write(G, coordinates, ss.str());
        }
#endif
    }
}

/**
 * Test whether points generated in hyperbolic space fulfill basic constraints
 */
TEST_F(HyperbolicGeneratorGTest, testHyperbolicPointGeneration) {
    count n = 1000;
    double stretch = Aux::Random::real(0.5, 1.5);
    double alpha = Aux::Random::real(0.5, 1.5);
    double R = HyperbolicSpace::hyperbolicAreaToRadius(n) * stretch;
    vector<double> angles(n, -1);
    vector<double> radii(n, -1);
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);
    for (index i = 0; i < n; i++) {
        EXPECT_GE(angles[i], 0);
        EXPECT_LT(angles[i], 2 * PI);
        EXPECT_GE(radii[i], 0);
        EXPECT_LE(radii[i], R);
    }
}

/**
 * Test whether the number edges generated by the hyperbolic generator agree at least roughly with
 * theory
 */
TEST_F(HyperbolicGeneratorGTest, testHyperbolicGenerator) {
    Aux::Random::setSeed(0, false);
    count n = 5000;
    double k = 16;
    count m = k * n / 2;
    HyperbolicGenerator gen(n, k, 7);
    Graph G = gen.generate();
    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_NEAR(G.numberOfEdges(), m, m / 5);
}

/**
 * Check consistency of graphs generated by the hyperbolic generator
 */
TEST_F(HyperbolicGeneratorGTest, testHyperbolicGeneratorConsistency) {
    Aux::Random::setSeed(0, false);
    count n = 5000;
    double k = 6;
    count m = n * k / 2;
    HyperbolicGenerator gen(n, k);
    Graph G = gen.generate();
    EXPECT_NEAR(G.numberOfEdges(), m, m / 5);
    ASSERT_TRUE(G.checkConsistency());
}

TEST_F(HyperbolicGeneratorGTest, testHyperbolicGeneratorMechanicGraphs) {
    Aux::Random::setSeed(0, false);
    count n = 2000;
    double k = 6;
    count m = n * k / 2;
    HyperbolicGenerator gen(n, k, 3, 0.14);
    Graph G = gen.generate();
    EXPECT_NEAR(G.numberOfEdges(), m, m / 10);
    ASSERT_TRUE(G.checkConsistency());
}

TEST_F(HyperbolicGeneratorGTest, debugHyperbolicHighTemperatureGraphs) {
    count n = 10000;
    double k = 10;
    double gamma = 3;
    count m = n * k / 2;
    for (double T = 0; T < 10; T += 0.1) {
        if (std::abs(T - 1) < 0.00001)
            continue;
        HyperbolicGenerator gen(n, k, gamma, T);
        Graph G = gen.generate();
        EXPECT_NEAR(G.numberOfEdges(), m, m / 10);
    }
}

TEST_F(HyperbolicGeneratorGTest, debugGiganticCollectionOfHyperbolicTemperatureGraphs) {
    for (index i = 0; i < 30; i++) {
        count n = 10000;
        double k = 10;
        double T = 0.1;
        count m = n * k / 2;
        HyperbolicGenerator gen(n, k, 3, T);
        Graph G = gen.generate();
        EXPECT_NEAR(G.numberOfEdges(), m, m / 10);
    }
}

TEST_F(HyperbolicGeneratorGTest, debugGiganticCollectionOfHyperbolicUnitDiskGraphs) {
    count n = 1000000;
    double k = 1;
    for (index i = 0; i < 7; i++) {
        count m = n * k / 2;
        HyperbolicGenerator gen(n, k, 7);
        Graph G = gen.generate();
        EXPECT_NEAR(G.numberOfEdges(), m, m / 5);
        EXPECT_TRUE(G.checkConsistency());
        k *= 2;
    }
}

} /* namespace NetworKit */

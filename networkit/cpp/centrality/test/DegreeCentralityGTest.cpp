// networkit-format
/*
 * DegreeCentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <array>

#include <gtest/gtest.h>

#include <networkit/centrality/DegreeCentrality.hpp>
#include <networkit/graph/Graph.hpp>

#include "CentralityTestFixture.hpp"

namespace NetworKit {

class DegreeCentralityGTest : public DirectedWeightedGTest {};

INSTANTIATE_TEST_SUITE_P(InstantiationName, DegreeCentralityGTest, directedWeightedParams());

TEST_P(DegreeCentralityGTest, testDegreeCentrality) {
    Graph g(8, false, isDirected());

    g.addEdge(0, 2);
    g.addEdge(0, 5);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 2);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(7, 7);

    DegreeCentrality dc(g, false, true, false);
    dc.run();
    EXPECT_EQ(dc.maximum(), 8);

    if (isDirected()) {
        std::array<int, 8> expectedResults{{2, 1, 3, 1, 1, 3, 0, 1}};
        for (long unsigned int i = 0; i < expectedResults.size(); i++)
            EXPECT_EQ(expectedResults[i], dc.score(i));
    } else {
        std::array<int, 8> expectedResults{{2, 1, 5, 2, 2, 6, 1, 2}};
        for (long unsigned int i = 0; i < expectedResults.size(); i++)
            EXPECT_EQ(expectedResults[i], dc.score(i));
    }
}

TEST_P(DegreeCentralityGTest, testInDegreeCentrality) {
    Graph g(8, false, isDirected());

    g.addEdge(0, 2);
    g.addEdge(0, 5);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 2);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(7, 7);

    DegreeCentrality dc(g, false, false, false);
    dc.run();

    if (isDirected()) {
        std::array<int, 8> expectedResults{{0, 0, 3, 1, 1, 4, 1, 2}};
        for (long unsigned int i = 0; i < expectedResults.size(); i++)
            EXPECT_EQ(expectedResults[i], dc.score(i));
    } else {
        std::array<int, 8> expectedResults{{2, 1, 5, 2, 2, 6, 1, 2}};
        for (long unsigned int i = 0; i < expectedResults.size(); i++)
            EXPECT_EQ(expectedResults[i], dc.score(i));
    }
}

TEST_P(DegreeCentralityGTest, testDegreeCentralityIgnoreSelfLoops) {
    Graph g(8, false, isDirected());

    g.addEdge(0, 2);
    g.addEdge(0, 5);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 2);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(7, 7);

    DegreeCentrality dc(g, false, true, true);
    dc.run();

    if (isDirected()) {
        std::array<int, 8> expectedResults{{2, 1, 2, 1, 1, 2, 0, 0}};
        for (long unsigned int i = 0; i < expectedResults.size(); i++)
            EXPECT_EQ(expectedResults[i], dc.score(i));
    } else {
        std::array<int, 8> expectedResults{{2, 1, 4, 2, 2, 5, 1, 1}};
        for (long unsigned int i = 0; i < expectedResults.size(); i++)
            EXPECT_EQ(expectedResults[i], dc.score(i));
    }
}

} /* namespace NetworKit */

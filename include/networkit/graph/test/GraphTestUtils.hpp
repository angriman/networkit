/*
 * GraphTestUtils.hpp
 *
 * Shared test helpers (custom gmock matchers) for the graph unit tests.
 */

#ifndef NETWORKIT_GRAPH_TEST_GRAPH_TEST_UTILS_HPP_
#define NETWORKIT_GRAPH_TEST_GRAPH_TEST_UTILS_HPP_

#include <gmock/gmock.h>

namespace NetworKit {

// Matches an edge-like value (anything with `.u` and `.v` members) against the expected endpoints.
MATCHER_P2(EdgeEq, expectedU, expectedV, "") {
    return arg.u == expectedU && arg.v == expectedV;
}

// Matches a weighted-edge-like value (`.u`, `.v`, `.weight`) against the expected endpoints and
// weight.
MATCHER_P3(WeightedEdgeEq, expectedU, expectedV, expectedWeight, "") {
    return arg.u == expectedU && arg.v == expectedV && arg.weight == expectedWeight;
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_TEST_GRAPH_TEST_UTILS_HPP_

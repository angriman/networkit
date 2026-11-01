// networkit-format
#ifndef NETWORKIT_CENTRALITY_TEST_CENTRALITY_TEST_FIXTURE_HPP_
#define NETWORKIT_CENTRALITY_TEST_CENTRALITY_TEST_FIXTURE_HPP_

#include <utility>

#include <gtest/gtest.h>

namespace NetworKit {

/**
 * Base fixture for tests that run on every combination of directed and weighted graphs.
 * The first parameter is weighted, the second one is directed.
 */
class DirectedWeightedGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    bool isWeighted() const noexcept { return GetParam().first; }
    bool isDirected() const noexcept { return GetParam().second; }
};

inline auto directedWeightedParams() {
    return testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                           std::make_pair(false, true), std::make_pair(true, true));
}

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_TEST_CENTRALITY_TEST_FIXTURE_HPP_

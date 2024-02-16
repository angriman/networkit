/*
 * OverlappingCommunityDetectionAlgorithm.cpp
 *
 *  Created on: 14.12.2020
 *      Author: John Gelhausen
 */

#include <stdexcept>
#include <networkit/base/Algorithm.hpp>
#include <networkit/community/OverlappingCommunityDetectionAlgorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

OverlappingCommunityDetectionAlgorithm::OverlappingCommunityDetectionAlgorithm(const Graph &G)
    : Algorithm(), G(&G), result(0) {
    // currently our community detection methods are not defined on directed graphs
    if (G.isDirected()) {
        throw std::runtime_error("This community detection method is undefined on directed graphs");
    }
}

const Cover &OverlappingCommunityDetectionAlgorithm::getCover() const {
    assureFinished();
    return result;
}

} /* namespace NetworKit */

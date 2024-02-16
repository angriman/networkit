/*
 * DissimilarityMeasure.cpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt
 */

#include <stdexcept>
#include <networkit/community/DissimilarityMeasure.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

double DissimilarityMeasure::getDissimilarity(const Graph &, const Cover &, const Cover &) {
    throw std::runtime_error("Not implemented");
}

} /* namespace NetworKit */

/*
 * EdgeCut.cpp
 *
 *  Created on: Jun 20, 2013
 *      Author: Henning
 */

#include <networkit/Globals.hpp>
#include <networkit/community/EdgeCut.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

double EdgeCut::getQuality(const Partition &zeta, const Graph &G) {
    double cutWeight = 0.0;
    G.forEdges([&](node u, node v, edgeweight w) {
        if (zeta[u] != zeta[v]) {
            cutWeight += w;
        }
    });
    return cutWeight;
}

} /* namespace NetworKit */

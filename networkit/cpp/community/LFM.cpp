/*
 * LFM.cpp
 *
 *  Created on: 10.12.2020
 *      Author: John Gelhausen
 */

#include <set>
#include <utility>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/LFM.hpp>
#include <networkit/community/OverlappingCommunityDetectionAlgorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/scd/SelectiveCommunityDetector.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

LFM::LFM(const Graph &G, SelectiveCommunityDetector &scd)
    : OverlappingCommunityDetectionAlgorithm(G), scd(&scd) {}

void LFM::run() {
    Aux::SignalHandler handler;

    Cover zeta(G->upperNodeIdBound());
    index o = 0;
    zeta.setUpperBound(o);

    G->forNodesInRandomOrder([&](node u) {
        handler.assureRunning();
        if (!zeta.contains(u)) {
            std::set<node> community = scd->expandOneCommunity(u);
            o++;
            zeta.setUpperBound(o);

            handler.assureRunning();

            for (node n : community) {
                zeta.addToSubset(o - 1, n);
            }
        }
    });

    result = std::move(zeta);
    hasRun = true;
}

} /* namespace NetworKit */

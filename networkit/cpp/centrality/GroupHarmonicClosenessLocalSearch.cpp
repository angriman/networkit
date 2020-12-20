/*
 *  GroupHarmonicClosenessLocalSearch.cpp
 *
 *  Created on: 20.12.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include "GroupHarmonicClosenessLocalSearchImpl.hpp"

#include <networkit/centrality/GroupHarmonicClosenessLocalSearch.hpp>

namespace NetworKit {
GroupHarmonicClosenessLocalSearch::GroupHarmonicClosenessLocalSearch(const Graph &G, count k)
    : weighted(G.isWeighted()) {
    if (weighted)
        implWeighted = std::make_unique<GroupHarmonicClosenessLocalSearchDetails::
                                            GroupHarmonicClosenessLocalSearchImplWeighted>(G, k);
    else
        implUnweighted =
            std::make_unique<GroupHarmonicClosenessLocalSearchDetails::
                                 GroupHarmonicClosenessLocalSearchImplUnweighted>(G, k);
}

GroupHarmonicClosenessLocalSearch::~GroupHarmonicClosenessLocalSearch() = default;

void GroupHarmonicClosenessLocalSearch::run() {
    if (weighted)
        implWeighted->run(maxSwaps);
    else
        implUnweighted->run(maxSwaps);
    hasRun = true;
}

const std::vector<node> &GroupHarmonicClosenessLocalSearch::groupMaxHarmonicCloseness() const {
    assureFinished();
    return weighted ? implWeighted->groupMaxHarmonicCloseness()
                    : implUnweighted->groupMaxHarmonicCloseness();
}

} // namespace NetworKit

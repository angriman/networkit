/*
 * WeightedDegreeCentrality.cpp
 *
 *  Created on: 23.10.2018
 *      Author: Eugenio Angriman
 */

#include "WeightedDegreeCentrality.h"
#include "../auxiliary/Parallel.h"

namespace NetworKit {

WeightedDegreeCentrality::WeightedDegreeCentrality(
    const Graph &G, const std::vector<double> &nodeScores)
    : G(G), nodeScores(nodeScores), n(G.upperNodeIdBound()) {
	hasRun = false;
}

void WeightedDegreeCentrality::run() {
	scoreData.resize(n, 0.);

	G.parallelForNodes([&](node u) {
		G.forNeighborsOf(
		    u, [&](node v) { scoreData[u] += G.weight(u, v) * nodeScores[v]; });
	});
	buildRanking();
	hasRun = true;
}

void WeightedDegreeCentrality::buildRanking() {
	rankingVec.resize(n);
	G.parallelForNodes(
	    [&](node u) { rankingVec[u] = std::make_pair(u, scoreData[u]); });
	Aux::Parallel::sort(
	    rankingVec.begin(), rankingVec.end(),
	    [&](std::pair<node, double> x, std::pair<node, double> y) {
		    if (x.second == y.second) {
			    return x.first < y.first;
		    }
		    return x.second > y.second;
	    });
}

std::vector<std::pair<node, double>> WeightedDegreeCentrality::ranking() const {
	if (!hasRun) {
		throw std::runtime_error("Call run method first.");
	}
	return rankingVec;
}

} /* namespace NetworKit */

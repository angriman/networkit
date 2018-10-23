/*
 * DegreeCentrality.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "DegreeCentrality.h"

namespace NetworKit {

DegreeCentrality::DegreeCentrality(const Graph &G, bool normalized, bool outDeg,
                                   bool ignoreSelfLoops, const bool weighted)
    : Centrality(G, normalized), outDeg(outDeg),
      ignoreSelfLoops(ignoreSelfLoops), weighted(weighted) {}

void DegreeCentrality::run() {
	scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);

	if (G.isDirected() && !outDeg) {
		G.parallelForNodes([&](node u) {
			scoreData[u] = weighted ? G.weightedDegreeIn(u) : G.degreeIn(u);
		});
	} else {
		G.parallelForNodes([&](node u) {
			scoreData[u] = weighted ? G.weightedDegree(u) : G.degree(u);
		});
	}

	if (normalized) {
		count maxDeg = maximum();
		G.parallelForNodes([&](node u) { scoreData[u] = scoreData[u] / maxDeg; });
	}
	hasRun = true;
}

double DegreeCentrality::maximum() {
	return G.numberOfNodes() - ignoreSelfLoops;
}

} /* namespace NetworKit */

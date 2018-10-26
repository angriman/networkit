/*
 * CurrentFlowCloseness.cpp
 *
 *  Created on: 01.10.2018
 *      Author: Alexander van der Grinten (alexander.vandergrinten@gmail.com)
 *      Adapted from Closeness.h
 */

#include "CurrentFlowCloseness.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/PrioQueue.h"
#include "../components/ConnectedComponents.h"
#include "../distance/CommuteTimeDistance.h"

namespace NetworKit {

CurrentFlowCloseness::CurrentFlowCloseness(
    const Graph &G, const std::vector<double> &nodeWeights, bool normalized,
    bool checkConnectedness)
    : Centrality(G, normalized), nodeWeights(nodeWeights) {
	if (checkConnectedness) {
		ConnectedComponents compo(G);
		compo.run();
		if (compo.numberOfComponents() != 1)
			throw std::runtime_error(
			    "CurrentFlowCloseness is undefined on disconnected graphs");
	}
}

void CurrentFlowCloseness::run() {
	scoreData.clear();
	scoreData.resize(G.numberOfNodes(), 0.0);

	edgeweight infDist = std::numeric_limits<edgeweight>::max();

	//	G.parallelForNodes([&](node s) {
	//		CommuteTimeDistance ctd{G};
	//		scoreData[s] = 1 / ctd.runSingleSource(s);
	//
	//	});

	CommuteTimeDistance ctd(G);
	std::cout << "Running ctd\n";
	ctd.run();
	G.forNodes([&](node s) {
		double sum = 0.0;
		G.forNodes([&](node t) {
			const double dist = ctd.distance(s, t);
			if (dist != 0 && dist != infDist) {
				sum += nodeWeights[s] / dist;
			}
		});
		scoreData[s] = sum;
	});

	if (normalized) {
		G.forNodes(
		    [&](node u) { scoreData[u] = scoreData[u] * (G.numberOfNodes() - 1); });
	}

	hasRun = true;
}

} /* namespace NetworKit */

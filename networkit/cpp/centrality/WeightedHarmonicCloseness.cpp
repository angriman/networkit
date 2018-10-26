/*
 * WeightedHarmonicCloseness.cpp
 *
 * Created on: 24.02.2018
 * 		 Author: Eugenio Angriman
 */

#include <memory>

#include "../distance/BFS.h"
#include "../distance/Dijkstra.h"
#include "../distance/SSSP.h"
#include "WeightedHarmonicCloseness.h"

namespace NetworKit {

WeightedHarmonicCloseness::WeightedHarmonicCloseness(
    const Graph &G, const std::vector<double> &nodeWeights,
    const std::vector<bool> &inGroup, bool normalized)
    : Centrality(G, normalized), nodeWeights(nodeWeights), inGroup(inGroup),
      n(G.upperNodeIdBound()) {}

void WeightedHarmonicCloseness::run() {
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	scoreData.clear();
	scoreData.resize(G.upperNodeIdBound(), 0.0);

	G.parallelForNodes([&](node v) {
		if (inGroup[v]) {
			return;
		}
		std::unique_ptr<SSSP> sssp;
		if (G.isWeighted()) {
			sssp.reset(new Dijkstra(G, v, true, true));
		} else {
			sssp.reset(new BFS(G, v, true, true));
		}

		sssp->run();

		std::vector<edgeweight> distances = sssp->getDistances();

		double sum = 0;
		for (node u = 0; u < n; ++u) {
			edgeweight dist = distances[u];
			if (dist != infDist && dist != 0) {
				sum += nodeWeights[u] / dist;
			}
		}
		scoreData[v] = sum;
	});

	if (normalized) {
		G.forNodes([&](node w) { scoreData[w] /= (G.numberOfNodes() - 1); });
	}

	hasRun = true;
}

double WeightedHarmonicCloseness::maximum() {
	return normalized ? (G.numberOfNodes() - 1) : 1.f;
}
} // namespace NetworKit

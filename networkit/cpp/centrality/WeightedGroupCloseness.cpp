#include <omp.h>

#include "../auxiliary/Log.h"
#include "../centrality/HarmonicCloseness.h"
#include "../distance/Dijkstra.h"
#include "WeightedGroupCloseness.h"

namespace NetworKit {

WeightedGroupCloseness::WeightedGroupCloseness(const Graph &G, const count k)
    : G(G), k(k), n(G.upperNodeIdBound()) {
	if (!G.isWeighted()) {
		throw std::runtime_error("Error: the graph is not weighted.");
	}

	if (k == 0 || k >= n) {
		throw std::runtime_error("Error: k must be within 0 and n-1.");
	}

	hasRun = false;
}

void WeightedGroupCloseness::run() {
	group.reserve(k);
	std::vector<bool> inGroup(n, false);

	HarmonicCloseness hc(G);
	hc.run();
	node topC = hc.ranking()[0].first;
	inGroup[topC] = true;
	group.push_back(topC);

	hasRun = true;
}

} // namespace NetworKit

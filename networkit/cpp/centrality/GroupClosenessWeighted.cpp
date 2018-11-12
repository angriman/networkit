#include "GroupClosenessWeighted.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/PrioQueue.h"

namespace NetworKit {

GroupClosenessWeighted::GroupClosenessWeighted(const Graph &G, const count k)
    : G(G), k(k), n(G.upperNodeIdBound()) {
	if (k == 0 || k > n) {
		throw std::runtime_error("Error: k must be at least 1 and at most n.");
	}
}

void GroupClosenessWeighted::init() {
	inGroup.assign(n, false);
	prio.assign(n, 0.0);
	group.reserve(k);
}

void GroupClosenessWeighted::computeInitialBound() {
	G.parallelForNodes([&](const node u) {
		prio[u] = 0;
		double curPrio = 0.0;
		double curDist;
		if (!inGroup[u]) {
			G.forNeighborsOf(u, [&](const node v, const edgeweight w) {
				curDist = dist[v];
				if (w < curDist) {
					curPrio -= curDist == infDist ? w : curDist - w;
				}
			});
			prio[u] = curPrio;
		}
	});
}

void GroupClosenessWeighted::run() {
	init();
	// TODO to save memory during runtime, safely destroy this object before
	// starting the algorithm.
	WeightedTopCloseness wtc(G, 1, false, false, true);
	wtc.run();

	group.push_back(wtc.topkNodesList()[0]);
	inGroup[group.back()] = true;
	dist = wtc.getTopNodeDist();
	tmpDist = dist;
	computeInitialBound();
	Aux::PrioQueue<double, node> Q(prio);

	double best = 0;

	while (group.size() < k) {

		break;
	}

	hasRun = true;
}
} // namespace NetworKit

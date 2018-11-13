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
	tmpDist.assign(n, infDist);
}

void GroupClosenessWeighted::computeInitialBound(const count &reachableFromTop,
                                                 const double &topSumDist) {
	G.parallelForNodes([&](const node u) {
		if (inGroup[u]) {
			return;
		}
		double curPrio = 0.0;
		double curDist;
		count newReachable = reachableFromTop;
		double newSum = topSumDist;
		G.forNeighborsOf(u, [&](const node v, const edgeweight w) {
			curDist = dist[v];
			if (curDist == infDist) {
				++newReachable;
				newSum += w;
			} else if (w < curDist) {
				newSum -= curDist - w;
			}
		});
		prio[u] = newSum * (n - 1.0) / (newReachable - 1.0) / (newReachable - 1.0);
	});
}

void GroupClosenessWeighted::bfsCut(const node &s) {}

void GroupClosenessWeighted::eraseVisitedEdges(
    std::vector<index> &topVisitedEdges) {
	for (auto eid : topVisitedEdges) {
		visitedEdges[eid] = true;
	}

	auto toErase = sortedEdges.begin();
	while (visitedEdges[toErase->first]) {
		++toErase;
	}
	sortedEdges.erase(sortedEdges.begin(), toErase);
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
	sortedEdges = wtc.getSortedEdges();
	visitedEdges = wtc.getVisitedEdgesVector();
	eraseVisitedEdges(wtc.getVisitedEdges());

	computeInitialBound(wtc.getTopNodeReachable(), wtc.getTopSum());

	double best = 0;

	while (group.size() < k) {
		Aux::PrioQueue<double, node> Q(prio);

		break;
	}

	hasRun = true;
}
} // namespace NetworKit

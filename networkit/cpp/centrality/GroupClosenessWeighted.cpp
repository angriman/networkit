#include "GroupClosenessWeighted.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Parallel.h"
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
	prio.assign(n, infDist);
	group.reserve(k);
	tmpDist.assign(n, infDist);
}

void GroupClosenessWeighted::computeInitialBound(
    const std::vector<double> &sortedDist, const double &sumDist) {
	G.parallelForNodes([&](const node u) {
		if (inGroup[u]) {
			return;
		}
		double imprUpperBound = 0.0, minDistance;
		double overlap = 0.0;
		count nearerNeigh = 0;
		double curDist, minWeightNeigh = infDist;
		G.forNeighborsOf(u, [&](const node v, const edgeweight w) {
			minWeightNeigh = std::min(minWeight, w);
			if (inGroup[v]) {
				overlap += 1.0;
			}
		});
		G.forNeighborsOf(u, [&](const node v, const edgeweight w) {
			curDist = dist[v];
			if (w < curDist) {
				imprUpperBound += std::min(w, minWeightNeigh + minWeight);
				++nearerNeigh;
			}
		});

		if (nearerNeigh > 0) {
			minDistance = minWeightNeigh + minWeight;
			INFO("Node ", u, " has ", minDistance, " min distance");
			count firstGreater = std::upper_bound(sortedDist.begin() + 1,
			                                      sortedDist.end(), minDistance) -
			                     sortedDist.begin();
			if (rL[u] < n - firstGreater) {
				firstGreater = n - rL[u];
			}
			double nonImprovableSum =
			    firstGreater < n / 2
			        ? std::accumulate(sortedDist.begin() + 1,
			                          sortedDist.begin() + firstGreater, 0.0)
			        : sumDist - std::accumulate(sortedDist.begin() + firstGreater,
			                                    sortedDist.end(), 0.0);

			imprUpperBound +=
			    nonImprovableSum +
			    std::max(0.0, n - firstGreater - G.degreeOut(u) + overlap - 1.0) *
			        minDistance;
		} else {
			prio[u] = sumDist;
		}
	});
}

void GroupClosenessWeighted::bfsCut(const node &s) { tmpDist[s] = 0.0; }

void GroupClosenessWeighted::run() {
	init();
	// TODO to save memory during runtime, safely destroy this object before
	// starting the algorithm.
	WeightedTopCloseness wtc(G, 1, false, false, true);
	wtc.run();

	group.push_back(wtc.topkNodesList()[0]);
	inGroup[group.back()] = true;
	dist = wtc.getTopNodeDist();
	std::vector<double> sortedDist(dist);
	Aux::Parallel::sort(sortedDist.begin(), sortedDist.end());
	reachable = wtc.getTopNodeReachable();
	rL = wtc.getReachL();
	minWeight = wtc.getMinWeight();
	computeInitialBound(sortedDist, wtc.getTopSum());

	node s;
	while (group.size() < k) {
		Aux::PrioQueue<double, node> Q(prio);
		while (Q.size() > 0) {
			s = Q.extractMin().second;
			INFO("Bfscut from ", s);
		}
		break;
	}

	hasRun = true;
}
} // namespace NetworKit

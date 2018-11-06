#include <cmath>
#include <omp.h>
#include <queue>

#include "../auxiliary/PrioQueue.h"
#include "../components/StronglyConnectedComponents.h"
#include "WeightedTopCloseness.h"

namespace NetworKit {
WeightedTopCloseness::WeightedTopCloseness(const Graph &G, const count k,
                                           const bool firstHeu,
                                           const bool secondHeu)
    : G(G), k(k), firstHeu(firstHeu), secondHeu(secondHeu),
      n(G.upperNodeIdBound()), infDist(std::numeric_limits<double>::max()) {
	if (k == 0 || k > n) {
		throw std::runtime_error("Error: k must be at least 1 and at most n.");
	}
}

void WeightedTopCloseness::init() {
	topkNodes.reserve(k);
	farness.assign(n, infDist);
	reachL.assign(n, 0);
	reachU.assign(n, 0);
}

void WeightedTopCloseness::computeReachable() {
	StronglyConnectedComponents sccs(G);
	sccs.run();

	const count N = sccs.numberOfComponents();
	Graph sccGraph(N, false, true);
	std::vector<bool> isAdj(N, false);
	std::vector<bool> reachFromMaxSCC(N, false);
	std::vector<count> reachLSCC(N, 0);
	std::vector<count> reachUSCC(N, 0);
	std::vector<count> reachUWithoutMaxSCC(N, 0);
	std::vector<bool> reachFromMaxScc(N, false);
	std::vector<bool> reachesMaxSCC(N, false);
	std::vector<std::vector<node>> sccVec(N);

	for (auto elem : sccs.getPartition().subsetSizeMap()) {
		sccVec[elem.first - 1].reserve(elem.second);
	}

	for (node u = 0; u < n; ++u) {
		sccVec[sccs.componentOfNode(u) - 1].push_back(u);
	}

	count maxSizeSCC = 0;
	// TODO parallelize
	for (count V = 0; V < N; ++V) {
		for (node u : sccVec[V]) {
			G.forNeighborsOf(u, [&](node v) {
				count W = sccs.componentOfNode(v) - 1;
				if (W != V && !isAdj[W]) {
					isAdj[W] = true;
					sccGraph.addEdge(V, W);
				}
			});
		}
		sccGraph.forNeighborsOf(V, [&](node adjSCC) { isAdj[adjSCC] = false; });
		if (sccGraph.degreeOut(V) > sccGraph.degreeOut(maxSizeSCC)) {
			maxSizeSCC = V;
		}
	}

	std::queue<node> Q;
	Q.push(maxSizeSCC);
	reachFromMaxSCC[maxSizeSCC] = true;
	count V;
	while (!Q.empty()) {
		V = Q.front();
		Q.pop();
		reachLSCC[maxSizeSCC] += sccVec[V].size();
		sccGraph.forNeighborsOf(V, [&](node W) {
			if (!reachFromMaxSCC[W]) {
				reachFromMaxSCC[W] = true;
				Q.push(W);
			}
		});
	}

	reachUSCC[maxSizeSCC] = reachLSCC[maxSizeSCC];
	reachesMaxSCC[maxSizeSCC] = true;

	for (node V = 0; V < N; ++V) {
		if (V == maxSizeSCC) {
			continue;
		}
		sccGraph.forNeighborsOf(V, [&](node W) {
			reachLSCC[V] = std::max(reachLSCC[V], reachLSCC[W]);
			if (!reachFromMaxSCC[W]) {
				reachUWithoutMaxSCC[V] += reachUWithoutMaxSCC[W];
			}
			reachUSCC[V] += reachUSCC[W];
			reachUSCC[V] = std::min(reachUSCC[V], n);
			reachesMaxSCC[V] = reachesMaxSCC[V] || reachesMaxSCC[W];
		});
		if (reachesMaxSCC[V]) {
			reachUSCC[V] += reachUWithoutMaxSCC[V];
		}

		reachLSCC[V] += sccVec[V].size();
		reachUSCC[V] += sccVec[V].size();
		reachUSCC[V] = std::min(reachUSCC[V], n);
	}

	for (node v = 0; v < n; ++v) {
		reachL[v] = reachLSCC[sccs.componentOfNode(v) - 1];
		reachU[v] = reachUSCC[sccs.componentOfNode(v) - 1];
	}
}

void WeightedTopCloseness::computeBounds() {}

void WeightedTopCloseness::run() {
	init();
	computeReachable();
	if (firstHeu) {
		DEBUG("Computing neighborhood-based lower bounds.");
		computeBounds();
	} else {
		DEBUG("Using weighted out degrees as lower bounds.");
		G.parallelForNodes([&](node u) {
			if (G.degreeOut(u) > 0) {
				farness[u] = -G.weightedDegree(u);
			}
		});
	}

	Aux::PrioQueue<double, node> q(farness);
	DEBUG("Done filling the queue.");
	hasRun = true;
}

} // namespace NetworKit

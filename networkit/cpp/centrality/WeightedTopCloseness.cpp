#include <cmath>
#include <omp.h>
#include <queue>

#include "../auxiliary/Parallel.h"
#include "../components/StronglyConnectedComponents.h"
#include "WeightedTopCloseness.h"

namespace NetworKit {
WeightedTopCloseness::WeightedTopCloseness(const Graph &G, const count k,
                                           const bool firstHeu,
                                           const bool secondHeu,
                                           const bool storeTopDist)
    : G(G), k(k), firstHeu(firstHeu), secondHeu(secondHeu),
      n(G.upperNodeIdBound()), storeTopDist(storeTopDist),
      minWeight(computeMinWeight()), kth(infDist), top(n) {
	if (k == 0 || k > n) {
		throw std::runtime_error("Error: k must be at least 1 and at most n.");
	}
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("Call index edges before.");
	}
}

void WeightedTopCloseness::init() {
	reachableFromTop = 0;
	topSum = 0.0;
	topkNodes.resize(k);
	farness.assign(n, infDist);
	reachL.assign(n, 0);
	dist.assign(n, infDist);
	if (storeTopDist) {
		topDist.assign(n, infDist);
	}
	reached.assign(n, false);
	lowerBoundDist.assign(n, infDist);
	nodesToReset.resize(n);
}

double WeightedTopCloseness::computeMinWeight() {
	double result = infDist;
	G.forEdges([&](const node u, const node v, const edgeweight w) {
		if (w < result) {
			result = w;
		}
	});
	return result;
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
	}
}

void WeightedTopCloseness::computeBounds() {
	// TODO: compute an alternative bound
	G.parallelForNodes([&](const node u) {
		if (G.degreeOut(u) > 0) {
			const double rL = reachL[u];
			double minWeightNeighbor = infDist;
			farness[u] = 0.0;
			G.forNeighborsOf(u, [&](const node v, const edgeweight w) {
				if (w < minWeightNeighbor) {
					minWeightNeighbor = w;
				}
			});
			G.forNeighborsOf(u, [&](const node v, const edgeweight w) {
				farness[u] += std::min(w, minWeightNeighbor + minWeight);
			});

			farness[u] += std::max(0.0, rL - 1.0 - G.degreeOut(u)) *
			              (minWeightNeighbor + minWeight);
			farness[u] *= (n - 1) / (rL - 1.0) / (rL - 1.0);
		} else {
			farness[u] = infDist;
		}
	});
}

void WeightedTopCloseness::bfsBound(const node &s) {}

double WeightedTopCloseness::bfsCut(const node &s) {
	reachedNodes = 1;
	nodesToReset[0] = s;
	reached[s] = true;
	dist[s] = 0.0;
	lowerBoundDist[s] = 0.0;
	// TODO use a heap and AVOID to reallocate at each bfscut!
	Aux::PrioQueue<double, node> pq(dist);
	const double &rL = reachL[s];
	d = 0.0;
	double sumDistLower = 0.0, newDist, lower = 0.0, distCur;
	node cur;
	std::pair<double, node> curPair;
	INFO("BFscut from ", s);

	while (pq.size() > 0) {
		curPair = pq.extractMin();
		cur = curPair.second;
		distCur = curPair.first;
		if (distCur == infDist) {
			// All reachable nodes have been visited, exiting the loop.
			break;
		}
		d += distCur;
		assert(lowerBoundDist[cur] != infDist);
		sumDistLower -= lowerBoundDist[cur];

		G.forNeighborsOf(cur, [&](const node v, const edgeweight w) {
			newDist = distCur + w;
			if (dist[v] > newDist) {
				dist[v] = newDist;
				pq.changeKey(newDist, v);
				if (!reached[v]) {
					reached[v] = true;
					nodesToReset[nodesToResetCount++] = v;
					++reachedNodes;
				} else {
					sumDistLower -= lowerBoundDist[v];
				}
				lowerBoundDist[v] = std::min(newDist, distCur + minWeight);
				sumDistLower += lowerBoundDist[v];
			}
		});

		lower = d + sumDistLower;
		if (pq.size() > 0) {
			// Next node to be picked, this distance is correct.
			curPair = pq.peekMin(0);
			newDist = curPair.first;
			cur = curPair.second;
			if (newDist != infDist) {
				lower += newDist - lowerBoundDist[cur];
			}
			if (rL > reachedNodes) {
				// FIXME what if newDist is inf??
				lower += (rL - reachedNodes) * (newDist + minWeight);
				lower *= (n - 1) / (rL - 1.0) / (rL - 1.0);
			} else {
				lower *= (n - 1) / (reachedNodes - 1.0) / (reachedNodes - 1.0);
			}
		} else {
			lower *= (n - 1) / (reachedNodes - 1.0) / (reachedNodes - 1.0);
		}

		if (lower >= kth) {
			break;
		}
	}

	return lower;
}

void WeightedTopCloseness::run() {
	init();
	computeReachable();
	computeBounds();

	node tmpNode;

	Aux::PrioQueue<double, node> Q(farness);
	DEBUG("Done filling the queue.");

	std::pair<double, node> topPair, bottomPair;
	node s;
	// TODO parallelize
	while (Q.size() > 0) {
		DEBUG("# of nodes to be analyzed: ", Q.size());
		auto p = Q.extractMin();
		s = p.second;

		INFO("Priority of node ", s, " is ", p.first);
		if (G.degreeOut(s) == 0 || farness[s] > kth) {
			INFO("Discarding node ", s);
			break;
		}

		if (secondHeu) {
			bfsBound(s);
		} else {
			nodesToResetCount = 1;
			farness[s] = bfsCut(s);
			INFO("Farness of ", s, " is ", farness[s]);
		}

		if (farness[s] < kth) {
			top.insert(-farness[s], s);
			// Pair with highest farness
			topPair = top.peekMin(0);

			if (storeTopDist) {
				// Node with lowest farness (i.e. highest closeness)
				tmpNode = top.peekMin(top.size() - 1).second;
				if (tmpNode == s) {
					// The current node is the one with highest closeness
					reachableFromTop = reachedNodes;
					topSum = d;
					std::copy(dist.begin(), dist.end(), topDist.begin());
				}
			}

			if (top.size() >= k) {
				while (top.size() > k) {
					top.extractMin();
				}
				kth = -topPair.first;
				DEBUG("new kth = ", kth);
			}
		}

		for (count i = 0; i < nodesToResetCount; ++i) {
			s = nodesToReset[i];
			dist[s] = infDist;
			lowerBoundDist[s] = infDist;
			reached[s] = false;
		}
	}

	// Filling result vector
	while (top.size() > 0) {
		topPair = top.extractMin();
		topkNodes[top.size()] = topPair.second;
	}

	hasRun = true;
}
} // namespace NetworKit

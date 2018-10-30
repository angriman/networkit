#include <cmath>

#include "../centrality/WeightedHarmonicCloseness.h"
#include "../distance/APSP.h"
#include "../distance/CommuteTimeDistance.h"
#include "IndependentCascade.h"
#include "LinearThreshold.h"
#include "TopCentrality.h"

namespace NetworKit {
TopCentrality::TopCentrality(const Graph &G, const count k, const bool useTopk)
    : G(G), k(k), n(G.upperNodeIdBound()), diam(computeDiamter()),
      useTopk(useTopk), hasRun(false), hasInitThreshold(false) {
	if (!G.isDirected() || !G.isWeighted()) {
		throw std::runtime_error("Error: the graph must be directed and weighted");
	}
}

void TopCentrality::run() {
	influencers.clear();
	influencers.reserve(k);
	inGroup.assign(n, false);
	nodeWeights.resize(n, 1.0);
	bool wasReversed = false;
	std::cout << "Running IC\n";
	Graph reversed = reverseWeights();
	WeightedHarmonicCloseness c(reversed, nodeWeights, inGroup, false);

	while (influencers.size() < k) {
		std::cout << "Algorithm: running " << influencers.size() << "/" << k
		          << std::endl;

		c.run();
		if (useTopk) {
			std::cout << "Using topk\n";
			for (auto elem : c.ranking()) {
				influencers.push_back(elem.first);
				if (influencers.size() == k) {
					break;
				}
			}
		} else {
			influencers.push_back(c.ranking()[0].first);
		}
		inGroup[influencers.back()] = true;
		nodeWeights[influencers.back()] = 0.;
		//		if (influencers.size() < k) {
		// This takes O(time_steps * n) time, we may want to not run it for
		// nothing.
		updateWeights();
		//		}
	}
	std::cout << "Algorithm has run\n";

	hasRun = true;
}

void TopCentrality::updateWeights() {
	bool stop = false;
	G.parallelForNodes([&](node u) { nodeWeights[u] = inGroup[u] ? 0.0 : 1.0; });
	// Prev and next are the probabilities that a node is influenced
	std::vector<double> prev(inGroup.begin(), inGroup.end());
	std::vector<double> next(n, 0.);

	// Prune the exploration once the probaibility of all nodes to be influenced
	// is below the threshold
	const double threshold = 1e-7;

	count iters = 0;
	while (!stop && iters++ <= 1.5 * diam) {
		++iters;
		stop = true;
		for (node w = 0; w < n; ++w) {
			if (inGroup[w]) {
				continue;
			}
			double pNotInfluencedNow = 1.;
			G.forInNeighborsOf(w, [&](node z) {
				pNotInfluencedNow *= 1. - G.weight(z, w) * prev[z];
			});
			if (pNotInfluencedNow < 1. - threshold) {
				stop = false;
			}
			double pInfluencedNow = (1. - pNotInfluencedNow) * nodeWeights[w];
			// Use log or a different function?
			nodeWeights[w] *= 1. - pInfluencedNow;
			next[w] = pInfluencedNow;
		}
		// Use references instead of copies
		std::copy(next.begin(), next.end(), prev.begin());
	}
	double sum = 0.0;
	for (double p : nodeWeights) {
		sum += 1.0 - p;
	}
	estimate.push_back(sum);
}

edgeweight TopCentrality::computeDiamter() const {
	Graph gUnw = G.toUnweighted();
	APSP apsp(gUnw);
	apsp.run();
	auto distances = apsp.getDistances();
	count maxDist = 0;
	const edgeweight infDist = std::numeric_limits<edgeweight>::max();

	for (auto curDist : distances) {
		for (edgeweight d : curDist) {
			if (d != infDist) {
				maxDist = d > maxDist ? d : maxDist;
			}
		}
	}
	return maxDist;
}

Graph TopCentrality::reverseWeights(const bool reverseDirection) const {
	Graph result(n, true, true);
	G.forEdges([&](node u, node v, edgeweight eid) {
		node from = reverseDirection ? v : u;
		node to = reverseDirection ? u : v;
		result.addEdge(from, to, 1.0 / G.weight(u, v));
	});
	return result;
}

Graph TopCentrality::getUndirected(const node s) const {
	Graph result(n, true, false);
	std::queue<node> q;
	q.push(s);
	std::vector<bool> visited(n, false);
	visited[s] = true;
	while (!q.empty()) {
		node u = q.front();
		q.pop();
		G.forNeighborsOf(u, [&](node v) {

		});
	}
	return result;
}

void TopCentrality::reverseEdges(Graph &graph, bool &wasReversed) {
	std::vector<std::pair<node, node>> edges(graph.numberOfEdges());
	std::vector<double> weights(graph.numberOfEdges(), 0);
	count idx = 0;
	graph.forEdges([&](node u, node v, edgeweight eid) {
		edges[idx] = std::make_pair(v, u);
		weights[idx++] = graph.weight(u, v);
	});
	idx = 0;
	for (auto e : edges) {
		graph.removeEdge(e.second, e.first);
		graph.addEdge(e.first, e.second, weights[idx++]);
	}
}

} // namespace NetworKit

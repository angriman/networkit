#include <cmath>

#include "../centrality/CurrentFlowCloseness.h"
#include "../centrality/WeightedHarmonicCloseness.h"
#include "../distance/APSP.h"
#include "IndependentCascade.h"
#include "TopCentrality.h"

namespace NetworKit {
TopCentrality::TopCentrality(const Graph &G, const count k, const Metric)
    : G(G), k(k), n(G.upperNodeIdBound()), diam(computeDiamter()), algo(algo),
      hasRun(false) {
	if (!G.isDirected() || !G.isWeighted()) {
		throw std::runtime_error("Error: the graph must be directed and weighted");
	}
}

void TopCentrality::run() {
	influencers.clear();
	influencers.reserve(k);
	inGroup.assign(n, false);
	nodeWeights.resize(n);
	Graph reversed = reverseWeights();
	while (influencers.size() < k) {
		if (algo == 0) {
			c.reset(
			    new WeightedHarmonicCloseness(reversed, nodeWeights, inGroup, false));
		} else if (algo == 1) {
			c.reset(new CurrentFlowCloseness(G, nodeWeights));
		} else {
			std::stringstream sstream;
			sstream << "Error, unsupported centrality metric: ";
			sstream << algo;
			throw std::runtime_error(sstream.str());
		}
		c->run();
		//		CurrentFlowCloseness c(G, nodeWeights);
		//		c.run();
		influencers.push_back(c->ranking()[0].first);
		inGroup[influencers.back()] = true;
		nodeWeights[influencers.back()] = 0.;
		if (influencers.size() < k) {
			// This takes O(time_steps * n) time, we may want to not run it for
			// nothing.
			updateWeights();
		}
	}

	hasRun = true;
}

void TopCentrality::updateWeights() {
	bool stop = false;
	G.parallelForNodes([&](node u) { nodeWeights[u] = inGroup[u] ? 0. : 1.; });
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
			nodeWeights[w] *= 1. - pInfluencedNow;
			next[w] = pInfluencedNow;
		}
		// Use references instead of copies
		std::copy(next.begin(), next.end(), prev.begin());
	}
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

Graph TopCentrality::reverseWeights() const {
	Graph result(n, true, true);
	G.forEdges([&](node u, node v, edgeweight eid) {
		result.addEdge(u, v, 1.0 / G.weight(u, v));
	});
	return result;
}

} // namespace NetworKit

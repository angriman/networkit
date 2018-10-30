#include <queue>

#include "../auxiliary/Parallel.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Random.h"
#include "../centrality/DegreeCentrality.h"
#include "../distance/APSP.h"
#include "IndependentCascade.h"
#include "InfluenceMaximization.h"

namespace NetworKit {
InfluenceMaximization::InfluenceMaximization(const Graph &G, const count k)
    : G(G), n(G.upperNodeIdBound()), k(k) {
	if (G.numberOfSelfLoops() > 0) {
		throw std::runtime_error("Error influence maximization does not support "
		                         "graphs with self loops.");
	}

	if (!G.isWeighted()) {
		throw std::runtime_error("Error, the graph is not weighted.");
	}

	hasRun = false;
}

void InfluenceMaximization::run() {
	computeDiameter();
	Aux::PrioQueue<double, node> top(n);
	nodeScores.assign(n, 1.);
	topInfuencers.clear();
	topInfuencers.reserve(k);

	DegreeCentrality dc(G);
	dc.run();
	auto ranking = dc.ranking();
	std::vector<double> bottom(n);
	bool stop = false;

	// TODO parallelize this loop
	while (!stop) {
		std::pair<node, double> r = ranking.front();
		ranking.erase(ranking.begin());

		// Estimating number of influenced nodes
		auto estimate = computeInflProb(r.first);
		double r_score = computeScore(estimate);
		if (top.size() < k) {
			if (top.size() == 0 || r_score < top.peekMin(0).first) {
				bottom = estimate;
			}
			top.insert(r_score, r.first);
		} else {
			stop = true;
		}
		updateNodeScores(estimate);
		// Can we only sort part of the vector?
		Aux::Parallel::sort(
		    ranking.begin(), ranking.end(),
		    [&](std::pair<node, double> x, std::pair<node, double> y) {
			    if (x.second == y.second) {
				    return x.first < y.first;
			    }
			    return x.second > y.second;
		    });
	}

	top.forElements(
	    [&](double key, node elem) { topInfuencers.push_back(elem); });
	hasRun = true;
}

void InfluenceMaximization::updateNodeScores(const std::vector<double> &probs) {
#pragma omp parallel for
	for (count i = 0; i < n; ++i) {
		nodeScores[i] *= 1. - probs[i];
	}
}

double
InfluenceMaximization::computeScore(const std::vector<double> &probs) const {
	double result = 0.;
#pragma omp parallel for reduction(+ : result)
	for (count i = 0; i < n; ++i) {
		result += probs[i];
	}
	return result;
}

std::vector<double> InfluenceMaximization::computeInflProb(const node &v) {
	bool stop = false;
	std::vector<double> prev(n, 0.);
	std::vector<double> next(n, 0.);
	std::vector<double> notInflAtPrev(n, 1.);
	std::vector<double> result(n, 1.);
	prev[v] = 1.;
	const double threshold = 1e-6;

	count iters = 0;
	while (!stop && iters++ <= 1.5 * diam) {
		stop = true;
		for (node w = 0; w < n; ++w) {
			if (w == v) {
				continue;
			}
			double pNotInfluencedNow = 1.;
			G.forInNeighborsOf(w, [&](node z) {
				pNotInfluencedNow *= 1. - G.weight(z, w) * prev[z];
			});
			if (pNotInfluencedNow < 1. - threshold) {
				stop = false;
			}
			double pInfluencedNow = (1. - pNotInfluencedNow) * notInflAtPrev[w];
			notInflAtPrev[w] *= 1. - pInfluencedNow;
			next[w] = pInfluencedNow;
			result[w] *= pNotInfluencedNow;
		}
		std::copy(next.begin(), next.end(), prev.begin());
	}
	for (count i = 0; i < n; ++i) {
		if (i != v) {
			result[i] = 1. - result[i];
		}
	}
	return result;
}

// TODO implement faster algorithm to approximate the diameter of directed
// graphs
count InfluenceMaximization::computeDiameter() {
	Graph gUnw = G.toUnweighted();
	APSP apsp(gUnw);
	apsp.run();
	auto distances = apsp.getDistances();
	edgeweight maxDist = 0;
	for (auto curDistances : distances) {
		std::sort(curDistances.begin(), curDistances.end(),
		          std::greater<edgeweight>());
		count i = 0;
		while (curDistances[i] > gUnw.numberOfNodes() && i < curDistances.size()) {
			++i;
		}
		if (i < curDistances.size()) {
			maxDist = std::max(maxDist, curDistances[i]);
		}
	}
	diam = maxDist;
	INFO("Diameter = ", diam);
}

} // namespace NetworKit

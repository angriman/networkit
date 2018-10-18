#include <queue>

#include "InfluenceMaximization.h"

namespace NetworKit {
InfluenceMaximization::InfluenceMaximization(const Graph &G)
    : G(G), n(G.upperNodeIdBound()) {}

std::vector<double> InfluenceMaximization::computeInflProb(const node &v) {
	bool stop = false;
	std::vector<double> prev(n, 0.);
	std::vector<double> next(n, 0.);
	std::vector<double> notInflAtPrev(n, 1.);
	std::vector<double> result(n, 1.);
	prev[v] = 1.;

	while (!stop) {
		stop = true;
		for (node w = 0; w < n; ++w) {
			double pNotInfluencedNow = 1.;
			G.forInNeighborsOf(w, [&](node z) {
				pNotInfluencedNow *= 1. - G.weight(z, w) * prev[z];
			});
			if (pNotInfluencedNow < 1.) {
				stop = false;
			}
			double pInfluencedNow = (1. - pNotInfluencedNow) * notInflAtPrev[w];
			notInflAtPrev[w] *= 1. - pInfluencedNow;
			next[w] = pInfluencedNow;
			result[w] *= pNotInfluencedNow;
		}
		INFO(next);
		std::copy(next.begin(), next.end(), prev.begin());
	}
	for (count i = 0; i < n; ++i) {
		result[i] = 1. - result[i];
	}

	return result;
}

std::vector<double>
InfluenceMaximization::performSimulation(const node &v, const count nIter) {
	const count max_threads = omp_get_max_threads();
	std::vector<std::vector<count>> influencedNodes(max_threads,
	                                                std::vector<count>(n, 0));
	std::vector<std::vector<bool>> influenced(max_threads,
	                                          std::vector<bool>(n, false));

#pragma omp parallel for
	for (omp_index i = 0; i < nIter; ++i) {
		std::vector<count> &nodes = influencedNodes[omp_get_thread_num()];
		std::vector<bool> &curInfl = influenced[omp_get_thread_num()];
		std::fill(curInfl.begin(), curInfl.end(), false);
		curInfl[v] = true;
		std::queue<node> activeNodes;
		activeNodes.push(v);
		while (!activeNodes.empty()) {
			node u = activeNodes.front();
			activeNodes.pop();
			G.forNeighborsOf(u, [&](node v) {
				if (!curInfl[v]) {
					if (Aux::Random::real(0, 1) <= G.weight(u, v)) {
						curInfl[v] = true;
						activeNodes.push(v);
					}
				}
			});
		}

		for (node u = 0; u < n; ++u) {
			if (curInfl[u]) {
				++nodes[u];
			}
		}
	}

	for (count i = 1; i < omp_get_max_threads(); ++i) {
		for (node u = 0; u < n; ++u) {
			influencedNodes[0][u] += influencedNodes[i][u];
		}
	}

	std::vector<double> result(n, 0.);
	for (node u = 0; u < n; ++u) {
		result[u] = (double)influencedNodes[0][u] / (double)nIter;
	}

	return result;
}
} // namespace NetworKit

#include <cmath>
#include <omp.h>

#include "../auxiliary/Random.h"
#include "IndependentCascade.h"

namespace NetworKit {

IndependentCascade::IndependentCascade(const Graph &G)
    : G(G), n(G.upperNodeIdBound()), max(0), min(G.upperNodeIdBound()),
      avg(0.0) {
	if (!G.isDirected() || !G.isWeighted()) {
		throw std::runtime_error(
		    "Error, the input graph is not directed and weighted.");
	}
}

void IndependentCascade::performSimulation(const std::vector<node> &seeds,
                                           const count nIter,
                                           const count randomSeed) {
	Aux::Random::setSeed(randomSeed, true);
	const count max_threads = omp_get_max_threads();
	std::vector<std::vector<count>> influencedNodes(max_threads,
	                                                std::vector<count>(n, 0));
	std::vector<std::vector<bool>> influenced(max_threads,
	                                          std::vector<bool>(n, false));
#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(nIter); ++i) {

		std::vector<count> &nodes = influencedNodes[omp_get_thread_num()];
		std::vector<bool> &curInfl = influenced[omp_get_thread_num()];

		std::fill(curInfl.begin(), curInfl.end(), false);
		std::queue<node> activeNodes;
		for (auto curNode : seeds) {
			activeNodes.push(curNode);
			curInfl[curNode] = true;
		}

		node u;
		count totalInfluencedNodes = 0;
		while (!activeNodes.empty()) {
			u = activeNodes.front();
			activeNodes.pop();
			G.forNeighborsOf(u, [&](node v) {
				if (!curInfl[v]) {
					if (Aux::Random::real(0, 1) <= G.weight(u, v)) {
						curInfl[v] = true;
						++nodes[v];
						activeNodes.push(v);
						++totalInfluencedNodes;
					}
				}
			});
		}

#pragma omp critical
		{
			max = totalInfluencedNodes > max ? totalInfluencedNodes : max;
			min = totalInfluencedNodes < min ? totalInfluencedNodes : min;
		}
	}

	for (count i = 1; i < max_threads; ++i) {
#pragma omp parallel for
		for (node u = 0; u < n; ++u) {
			influencedNodes[0][u] += influencedNodes[i][u];
		}
	}

	avg = 0.0;
#pragma omp parallel for reduction(+ : avg)
	for (node u = 0; u < n; ++u) {
		avg += (double)influencedNodes[0][u] / (double)nIter;
	}
	hasRun = true;
}
} // namespace NetworKit

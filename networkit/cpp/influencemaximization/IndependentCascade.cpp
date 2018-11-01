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
#pragma omp parallel for schedule(dynamic)
	for (omp_index i = 0; i < static_cast<omp_index>(nIter); ++i) {

		const count threadId = omp_get_thread_num();
		std::vector<count> &threadInflNodes = influencedNodes[threadId];
		std::vector<bool> &threadInfluenced = influenced[threadId];

		std::fill(threadInfluenced.begin(), threadInfluenced.end(), false);
		std::queue<node> activeNodes;
		for (auto curNode : seeds) {
			activeNodes.push(curNode);
			threadInfluenced[curNode] = true;
		}

		node u;
		count totalInfluencedNodes = seeds.size();
		while (!activeNodes.empty()) {
			u = activeNodes.front();
			activeNodes.pop();
			G.forNeighborsOf(u, [&](node v) {
				if (!threadInfluenced[v]) {
					if (Aux::Random::real(0, 1) <= G.weight(u, v)) {
						threadInfluenced[v] = true;
						++threadInflNodes[v];
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

	avg = 0.0;
	for (count i = 0; i < max_threads; ++i) {
		for (node u = 0; u < n; ++u) {
			avg += influencedNodes[i][u];
			if (std::find(seeds.begin(), seeds.end(), u) != seeds.end()) {
				assert(influencedNodes[i][u] == 0);
			}
		}
	}

	avg /= (double)nIter;
	avg += seeds.size();

	hasRun = true;
}
} // namespace NetworKit

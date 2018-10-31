#include <cmath>
#include <omp.h>

#include "../auxiliary/Log.h"
#include "../auxiliary/Parallel.h"
#include "../centrality/HarmonicCloseness.h"
#include "../distance/Dijkstra.h"
#include "WeightedGroupCloseness.h"

namespace NetworKit {

WeightedGroupCloseness::WeightedGroupCloseness(const Graph &G, const count k)
    : G(G), k(k), n(G.upperNodeIdBound()) {
	if (!G.isWeighted()) {
		throw std::runtime_error("Error: the graph is not weighted.");
	}

	if (k == 0 || k >= n) {
		throw std::runtime_error("Error: k must be within 0 and n-1.");
	}

	hasRun = false;
}

void WeightedGroupCloseness::run() {
	const edgeweight infDist = std::numeric_limits<edgeweight>::max();
	newDistances.resize(omp_get_max_threads(),
	                    std::vector<edgeweight>(n, infDist));
	distances.assign(n, infDist);
	group.reserve(k);
	inGroup.assign(n, false);
	margGain.assign(n, 0);

	while (group.size() < k) {
		updateMarginalGain();
	}
	hasRun = true;
}

void WeightedGroupCloseness::updateMarginalGain() {
	std::vector<std::pair<node, double>> ranking;
	ranking.resize(n);
	G.parallelForNodes([&](node s) {
		if (inGroup[s]) {
			ranking[s] = std::make_pair(s, 0);
			return;
		}

		Dijkstra dijkstra(G, s, false);
		dijkstra.run();

		std::vector<edgeweight> &curNewDistances =
		    newDistances[omp_get_thread_num()];
		curNewDistances = dijkstra.getDistances();
		double curMargGain = 0;
		for (node u = 0; u < n; ++u) {
			curNewDistances[u] = std::min(curNewDistances[u], distances[u]);
			if (curNewDistances[u] != 0 && !inGroup[u]) {
				curMargGain += 1.0 / curNewDistances[u];
			}
		}
		ranking[s] = std::make_pair(s, curMargGain);
	});

	Aux::Parallel::sort(
	    ranking.begin(), ranking.end(),
	    [&](std::pair<node, double> x, std::pair<node, double> y) {
		    if (x.second == y.second) {
			    return x.first < y.first;
		    }
		    return x.second > y.second;
	    });

	group.push_back(ranking.front().first);
	inGroup[group.back()] = true;
	Dijkstra dijkstra(G, group.back(), false);
	dijkstra.run();
	auto curNewDistances = dijkstra.getDistances();
	G.parallelForNodes([&](node s) {
		distances[s] = std::min(curNewDistances[s], distances[s]);
		if (inGroup[s]) {
			assert(distances[s] == 0);
		}
	});
}

} // namespace NetworKit

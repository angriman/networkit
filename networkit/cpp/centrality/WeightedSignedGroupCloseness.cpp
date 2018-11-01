#include <cmath>
#include <omp.h>

#include "../auxiliary/Log.h"
#include "../auxiliary/Parallel.h"
#include "../centrality/HarmonicCloseness.h"
#include "../distance/Dijkstra.h"
#include "WeightedSignedGroupCloseness.h"

namespace NetworKit {

WeightedSignedGroupCloseness::WeightedSignedGroupCloseness(const Graph &G,
                                                           const count k)
    : WeightedGroupCloseness(G, k) {
	if (!G.isWeighted()) {
		throw std::runtime_error("Error: the graph is not weighted.");
	}

	if (k == 0 || k >= n) {
		throw std::runtime_error("Error: k must be within 0 and n-1.");
	}

	hasRun = false;
}

void WeightedSignedGroupCloseness::run() {
	distances.assign(n, infDist);
	spSign.assign(n, 1);
	group.reserve(k);
	inGroup.assign(n, false);
	contributionSign.resize(omp_get_max_threads(), std::vector<int>(n, 1));

	while (group.size() < k) {
		updateMarginalGain();
	}

	hasRun = true;
}

void WeightedSignedGroupCloseness::updateMarginalGain() {
	auto signFunc = [&](const edgeweight w) -> int { return w >= 0 ? 1 : -1; };
	std::vector<std::pair<node, double>> ranking;
	ranking.reserve(n - group.size());

	G.parallelForNodes([&](node s) {
		if (inGroup[s]) {
			return;
		}

		Dijkstra dijkstra(G, s, true);
		dijkstra.run();
		auto curNewDistances = dijkstra.getDistances();

		std::vector<bool> processed(n, false);
		auto &sign = contributionSign[omp_get_thread_num()];
		sign[s] = 1;

		double curMargGain = 0;
		G.forNodes([&](const node t) {
			if (s == t || curNewDistances[t] >= distances[t] || inGroup[t]) {
				return;
			}
			auto path = dijkstra.getPath(t);
			if (path.size() > 0) {
				auto source = path.begin();
				auto target = path.begin() + 1;
				while (target != path.end()) {
					if (!processed[*target]) {
						sign[*target] =
						    signFunc(sign[*source] * G.weight(*source, *target));
						processed[*target] = true;
					}
					++source;
					++target;
				}
				curMargGain += sign[t] / curNewDistances[t];
			}
		});

#pragma omp critical
		{ ranking.push_back(std::make_pair(s, curMargGain)); }
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
	distances[group.back()] = 0;

	Dijkstra dijkstra(G, group.back(), true);
	dijkstra.run();
	auto curNewDistances = dijkstra.getDistances();

	G.parallelForNodes([&](node t) {
		distances[t] = std::min(distances[t], curNewDistances[t]);
	});
}
} // namespace NetworKit

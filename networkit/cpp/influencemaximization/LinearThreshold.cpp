#include <queue>
#include <vector>

#include "LinearThreshold.h"

namespace NetworKit {
LinearThreshold::LinearThreshold(const Graph &G,
                                 const std::vector<double> &threshold)
    : G(G), n(G.upperNodeIdBound()), threshold(threshold), hasRun(false) {
	if (!G.isWeighted()) {
		throw std::runtime_error("Error: the graph is not weighted.");
	}
}

void LinearThreshold::performSimulation(const std::vector<node> &seeds) {
	std::vector<bool> influenced(n, false);
	std::queue<node> queue;
	influencedNodes = 0;
	std::vector<double> curThreshold(n);
	std::copy(threshold.begin(), threshold.end(), curThreshold.begin());
	std::cout << "Size " << threshold.size() << " first = " << threshold[0]
	          << std::endl;
	for (auto curNode : seeds) {
		queue.push(curNode);
		influenced[curNode] = true;
	}

	node u;
	while (!queue.empty()) {
		u = queue.front();
		queue.pop();
		G.forNeighborsOf(u, [&](node v) {
			if (!influenced[v]) {
				curThreshold[v] -= (double)G.weight(u, v);
				if (curThreshold[v] <= 0) {
					influenced[v] = true;
					++influencedNodes;
					queue.push(v);
				}
			}
		});
	}

	hasRun = true;
}

} // namespace NetworKit

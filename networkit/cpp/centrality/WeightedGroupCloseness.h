#ifndef WEIGHTEDGROUPCLOSENESS_H_
#define WEIGHTEDGROUPCLOSENESS_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

class WeightedGroupCloseness : public Algorithm {

public:
	WeightedGroupCloseness(const Graph &G, const count k = 1);
	void run() override;
	std::vector<node> groupMaxCloseness();

protected:
	const Graph &G;
	const count k;
	const count n;
	const edgeweight infDist = std::numeric_limits<edgeweight>::max();
	std::vector<node> group;
	std::vector<edgeweight> distances;
	std::vector<bool> inGroup;

	void updateMarginalGain();
	std::vector<edgeweight> computeDistances(const node &s);
};

inline std::vector<node> WeightedGroupCloseness::groupMaxCloseness() {
	assureFinished();
	return group;
}

} // namespace NetworKit

#endif

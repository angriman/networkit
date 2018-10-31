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
	std::vector<node> group;
};

inline std::vector<node> WeightedGroupCloseness::groupMaxCloseness() {
	assureFinished();
	return groupMaxCloseness();
}

} // namespace NetworKit

#endif

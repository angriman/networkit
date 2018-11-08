#ifndef GROUPCLOSENESSWEIGHTED_H_
#define GROUPCLOSENESSWEIGHTED_H_

#include "../base/Algorithm.h"
#include "../centrality/WeightedTopCloseness.h"
#include "../graph/Graph.h"

namespace NetworKit {

class GroupClosenessWeighted : public Algorithm {

public:
	GroupClosenessWeighted(const Graph &G, const count k = 1);
	void run() override;
	std::vector<node> groupMaxCloseness() const;

private:
	const Graph &G;
	const count k;
	const count n;
	const double infDist = std::numeric_limits<double>::max();
	std::vector<bool> inGroup;
	std::vector<node> group;
	std::vector<double> dist;
	std::vector<double> tmpDist;
	std::vector<double> prio;

	void init();
	void computeInitialBound();
};

inline std::vector<node> GroupClosenessWeighted::groupMaxCloseness() const {
	assureFinished();
	return group;
}
} // namespace NetworKit

#endif

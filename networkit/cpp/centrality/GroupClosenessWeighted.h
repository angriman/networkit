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
	count reachable;
	double minWeight;
	std::vector<bool> inGroup;
	std::vector<node> group;
	std::vector<double> dist;
	std::vector<double> tmpDist;
	std::vector<double> prio;
	std::vector<count> rL;

	void init();
	void computeInitialBound(const std::vector<double> &sortedDist,
	                         const double &sumDist);
	void bfsCut(const node &s);
};

inline std::vector<node> GroupClosenessWeighted::groupMaxCloseness() const {
	assureFinished();
	return group;
}
} // namespace NetworKit

#endif

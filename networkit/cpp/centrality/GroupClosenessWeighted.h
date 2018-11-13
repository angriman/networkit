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
	std::vector<std::pair<index, double>> sortedEdges;
	std::vector<bool> visitedEdges;

	void init();
	void computeInitialBound(const count &reachableFromTop,
	                         const double &topSumDist);
	void bfsCut(const node &s);
	void eraseVisitedEdges(std::vector<index> &topVisitedEdges);
};

inline std::vector<node> GroupClosenessWeighted::groupMaxCloseness() const {
	assureFinished();
	return group;
}
} // namespace NetworKit

#endif

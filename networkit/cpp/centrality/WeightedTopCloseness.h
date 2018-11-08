#ifndef WEIGHTEDTOPCLOSENESS_H_
#define WEIGHTEDTOPCLOSENESS_H_

#include "../auxiliary/Log.h"
#include "../auxiliary/PrioQueue.h"
#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {
class WeightedTopCloseness : public Algorithm {
public:
	WeightedTopCloseness(const Graph &G, const count k = 1,
	                     const bool firstHeu = true, const bool secondHeu = true,
	                     const bool storeTopDist = false);
	void run() override;
	std::vector<node> topkNodesList() const;
	std::vector<double> getTopNodeDist() const;

private:
	const Graph &G;
	const count k;
	const bool firstHeu;
	const bool secondHeu;
	const count n;
	const bool storeTopDist;
	const double infDist = std::numeric_limits<double>::max();
	double kth;
	node tmpNode, topNode;
	std::vector<node> topkNodes;
	std::vector<double> farness;
	std::vector<bool> reached;
	std::vector<count> reachL;
	std::vector<double> dist;
	std::vector<double> topDist;
	std::vector<double> lowerBoundDist;
	std::vector<std::pair<index, double>> sortedEdges;
	std::vector<bool> visitedEdges;
	std::vector<node> nodesToReset;
	std::vector<index> edgesToReset;

	Aux::PrioQueue<double, node> top;

	void init();
	void computeReachable();
	void computeBounds();
	void bfsBound(const node &s);
	double bfsCut(const node &s);
	double minUnvisitedEdge() const;
};

inline std::vector<node> WeightedTopCloseness::topkNodesList() const {
	assureFinished();
	return topkNodes;
}

inline std::vector<double> WeightedTopCloseness::getTopNodeDist() const {
	if (!storeTopDist) {
		throw std::runtime_error("Error, storeTopDist should be set to true");
	}
	return topDist;
}
} // namespace NetworKit

#endif

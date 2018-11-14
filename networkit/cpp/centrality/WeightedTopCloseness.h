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
	count getTopNodeReachable() const;
	double getTopSum() const;
	std::vector<count> &getReachL() {
		assureFinished();
		return reachL;
	}
	double getMinWeight() const;

private:
	const Graph &G;
	const count k;
	const bool firstHeu;
	const bool secondHeu;
	const count n;
	const bool storeTopDist;
	const double infDist = std::numeric_limits<double>::max();
	const double minWeight;
	double kth;
	double topSum;
	double d;
	count reachedNodes, reachableFromTop;
	count nodesToResetCount;

	std::vector<node> topkNodes;
	std::vector<double> farness;
	std::vector<bool> reached;
	std::vector<count> reachL;
	std::vector<double> dist;
	std::vector<double> topDist;
	std::vector<double> lowerBoundDist;
	std::vector<node> nodesToReset;

	Aux::PrioQueue<double, node> top;

	void init();
	void computeReachable();
	void computeBounds();
	void bfsBound(const node &s);
	double bfsCut(const node &s);
	bool checkStoreTopDist() const;
	double computeMinWeight();
};

inline std::vector<node> WeightedTopCloseness::topkNodesList() const {
	assureFinished();
	return topkNodes;
}

inline bool WeightedTopCloseness::checkStoreTopDist() const {
	if (!storeTopDist) {
		throw std::runtime_error("Error, storeTopDist should be set to true");
	}
}
inline std::vector<double> WeightedTopCloseness::getTopNodeDist() const {
	checkStoreTopDist();
	assureFinished();
	return topDist;
}

inline count WeightedTopCloseness::getTopNodeReachable() const {
	checkStoreTopDist();
	assureFinished();
	return reachableFromTop;
}

inline double WeightedTopCloseness::getTopSum() const {
	checkStoreTopDist();
	assureFinished();
	return topSum;
}

inline double WeightedTopCloseness::getMinWeight() const {
	assureFinished();
	return minWeight;
}
} // namespace NetworKit
#endif

#ifndef WEIGHTEDTOPCLOSENESS_H_
#define WEIGHTEDTOPCLOSENESS_H_

#include "../auxiliary/Log.h"
#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {
class WeightedTopCloseness : public Algorithm {
public:
	WeightedTopCloseness(const Graph &G, const count k = 1,
	                     const bool firstHeu = true, const bool secondHeu = true);
	void run() override;
	std::vector<node> topkNodesList() const;

private:
	const Graph &G;
	const count k;
	const bool firstHeu;
	const bool secondHeu;
	const count n;
	const double infDist;
	std::vector<node> topkNodes;
	std::vector<double> farness;
	std::vector<count> reachL;
	std::vector<count> reachU;

	void init();
	void computeReachable();
	void computeBounds();
};

inline std::vector<node> WeightedTopCloseness::topkNodesList() const {
	assureFinished();
	return topkNodes;
}
} // namespace NetworKit

#endif

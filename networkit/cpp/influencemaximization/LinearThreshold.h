#ifndef LINEAR_THRESHOLD_H_
#define LINEAR_THRESHOLD_H_

#include "../graph/Graph.h"

namespace NetworKit {

class LinearThreshold {
public:
	LinearThreshold(const Graph &G, const std::vector<double> &threshold);
	void performSimulation(const std::vector<node> &seeds);
	count getInfluencedNodes() const;

private:
	const Graph &G;
	const std::vector<double> &threshold;
	const count n;
	count influencedNodes;
	void checkHasRun() const;
	bool hasRun;
};

inline void LinearThreshold::checkHasRun() const {
	if (!hasRun) {
		throw std::runtime_error("Call performSimulation method first.");
	}
}

inline count LinearThreshold::getInfluencedNodes() const {
	checkHasRun();
	return influencedNodes;
}

} // namespace NetworKit

#endif

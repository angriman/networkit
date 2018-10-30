#ifndef TOPCENTRALITY_H_
#define TOPCENTRALITY_H_

#include "../base/Algorithm.h"
#include "../centrality/Centrality.h"
#include "../graph/Graph.h"

namespace NetworKit {

class TopCentrality : public Algorithm {
public:
	TopCentrality(const Graph &G, const count k = 1, const bool useTopk = false);
	void run() override;
	std::vector<node> getInfluencers() const;
	std::vector<double> getEstimate() const { return estimate; }
	void setThreshold(const std::vector<double> threshold) {
		if (threshold.size() != n) {
			throw std::runtime_error("Error: the length of the threshold vector must "
			                         "be equal to the number of nods.");
		}
		hasInitThreshold = true;
		this->threshold = threshold;
	};

private:
	const Graph &G;
	const count k;
	const count n;
	const edgeweight diam;
	const bool useTopk;
	std::vector<node> influencers;
	std::vector<double> nodeWeights;
	std::vector<bool> inGroup;
	void checkHasRun() const;
	void updateWeights();
	void runLT();
	edgeweight computeDiamter() const;
	Graph reverseWeights(const bool reverseDirection = false) const;
	std::vector<double> estimate;
	Graph getUndirected(const node s) const;
	void reverseEdges(Graph &graph, bool &wasReversed);
	bool hasRun, hasInitThreshold;
	std::vector<double> threshold;
};

inline void TopCentrality::checkHasRun() const {
	if (!hasRun) {
		throw std::runtime_error("Call run method first.");
	}
}

inline std::vector<node> TopCentrality::getInfluencers() const {
	checkHasRun();
	return influencers;
}
} // namespace NetworKit
#endif

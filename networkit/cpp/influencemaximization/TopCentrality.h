#ifndef TOPCENTRALITY_H_
#define TOPCENTRALITY_H_

#include "../base/Algorithm.h"
#include "../centrality/CurrentFlowCloseness.h"
#include "../centrality/WeightedHarmonicCloseness.h"
#include "../graph/Graph.h"

namespace NetworKit {
class TopCentrality : public Algorithm {
public:
	TopCentrality(const Graph &G, const count k = 1);
	void run() override;
	std::vector<node> getInfluencers() const;

private:
	const Graph &G;
	const count k;
	const count n;
	const edgeweight diam;
	std::vector<node> influencers;
	std::vector<double> nodeWeights;
	std::vector<bool> inGroup;
	void checkHasRun() const;
	void updateWeights();
	edgeweight computeDiamter() const;
	Graph reverseWeights() const;
	bool hasRun;
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

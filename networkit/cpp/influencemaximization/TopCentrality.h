#ifndef TOPCENTRALITY_H_
#define TOPCENTRALITY_H_

#include "../base/Algorithm.h"
#include "../centrality/Centrality.h"
#include "../graph/Graph.h"

namespace NetworKit {

enum Metric { CLOSENESS = 0, CURRENT_FLOW = 1, KATZ = 2 };
class TopCentrality : public Algorithm {
public:
	TopCentrality(const Graph &G, const count k = 1,
	              const Metric algo = CLOSENESS);
	void run() override;
	std::vector<node> getInfluencers() const;

private:
	const Graph &G;
	const count k;
	const count n;
	const edgeweight diam;
	const Metric algo;
	std::unique_ptr<Centrality> c;
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

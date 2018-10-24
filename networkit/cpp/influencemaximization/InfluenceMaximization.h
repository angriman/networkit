#ifndef INFLUENCE_MAXIMIZATION_H_
#define INFLUENCE_MAXIMIZATION_H_

#include <omp.h>
#include <vector>

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

class InfluenceMaximization : public Algorithm {
public:
	InfluenceMaximization(const Graph &G, const count k = 1);
	std::vector<double> computeInflProb(const node &v);
	void run();
	std::vector<node> getTopInfluencers() const;
	count computeDiameter();

private:
	const Graph &G;
	const count n, k;
	count diam;

	std::vector<double> nodeScores;
	std::vector<node> topInfuencers;

	double computeScore(const std::vector<double> &probs) const;
	void updateNodeScores(const std::vector<double> &probs);
};

inline std::vector<node> InfluenceMaximization::getTopInfluencers() const {
	if (!hasRun) {
		throw std::runtime_error("Call run method first");
	}
	return topInfuencers;
}
} // namespace NetworKit
#endif

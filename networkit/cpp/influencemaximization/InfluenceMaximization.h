#ifndef INFLUENCE_MAXIMIZATION_H_
#define INFLUENCE_MAXIMIZATION_H_

#include <omp.h>
#include <vector>

#include "../graph/Graph.h"

namespace NetworKit {

class InfluenceMaximization {
public:
	InfluenceMaximization(const Graph &G);
	std::vector<double> computeInflProb(const node &v);
	std::vector<double> performSimulation(const node &v, const count nIter = 1e6);
	count computeDiameter();

private:
	const Graph &G;
	const count n;
	count diam;
};
} // namespace NetworKit
#endif

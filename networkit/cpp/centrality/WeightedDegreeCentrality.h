/*
 * WeightedDegreeCentrality.h
 *
 *  Created on: 23.10.2018
 *      Author: Eugenio Angriman
 */

#ifndef WEIGHTEDDEGREECENTRALITY_H_
#define WEIGHTEDDEGREECENTRALITY_H_

#include "../graph/Graph.h"

namespace NetworKit {

class WeightedDegreeCentrality {
public:
	WeightedDegreeCentrality(const Graph &G,
	                         const std::vector<double> &nodeScores);

	void run();
	std::vector<std::pair<node, double>> ranking() const;

private:
	const Graph &G;
	const std::vector<double> &nodeScores;
	const count n;
	std::vector<double> scoreData;
	bool hasRun;

	void buildRanking();
	std::vector<std::pair<node, double>> rankingVec;
};
} /* namespace NetworKit */
#endif /* WEIGHTEDDEGREECENTRALITY_H_ */

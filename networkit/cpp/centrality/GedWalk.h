/* GroupCloseness.h
 *
 *  Created on: 18.6.2018
 *      Author: Eugenio Angriman
 */

#ifndef GEDWALK_H_
#define GEDWALK_H_

#include "../algebraic/CSRMatrix.h"
#include "../auxiliary/PrioQueue.h"
#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

class GedWalk : public Algorithm {

protected:
	Graph &G;
	const count n, k;
	const double alpha;
	double oldScore, newScore;
	const uint8_t maxLevel = 10;
	const std::vector<double> alphas;

	std::vector<node> group;
	std::vector<bool> inGroup;
	std::vector<bool> isExact;
	std::vector<std::vector<count>> pathsIn, pathsOut;
	std::vector<double> score;
	Aux::PrioQueue<double, node> q;
	CSRMatrix gMat;

	void init();
	void computePaths();
	void computeUpperBound();
	void computeMarginalGain(const node &z);
	double computeScore();
	void checkHasRun() const;
	count maxDegree() const;

public:
	GedWalk(Graph &G, const count k = 1, const double alpha = -1.);
	std::vector<node> groupMaxGedWalk();
	void run() override;
	double getApxScore() const;
	double scoreOfGroup(const std::vector<node> &set);
	std::vector<double> initAlphas() const;
};

inline void GedWalk::checkHasRun() const {
	if (!hasRun) {
		throw std::runtime_error("Run method has not been called");
	}
}

inline count GedWalk::maxDegree() const {
	count result = 0;
#pragma omp parallel for reduction(max : result)
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		count curDeg = G.degree(i);
		result = result < curDeg ? curDeg : result;
	}
	return result;
}

inline std::vector<node> GedWalk::groupMaxGedWalk() { return group; }

inline double GedWalk::scoreOfGroup(const std::vector<node> &set) {
	if (!hasRun) {
		init();
	}
	std::fill(pathsIn[0].begin(), pathsIn[0].end(), 0);
	std::fill(pathsOut[0].begin(), pathsOut[0].end(), 1);
	std::fill(inGroup.begin(), inGroup.end(), false);
	for (auto i = set.begin(); i < set.end(); ++i) {
		inGroup[*i] = true;
		pathsIn[0][*i] = 1;
		pathsOut[0][*i] = 0;
	}
	computePaths();
	return computeScore();
}

inline double GedWalk::getApxScore() const {
	checkHasRun();
	return newScore;
}

inline std::vector<double> GedWalk::initAlphas() const {
	std::vector<double> result(maxLevel - 1);
	for (uint8_t i = 0; i < result.size(); ++i) {
		result[i] = std::pow(alpha, i + 1);
	}

	return result;
}

} // namespace NetworKit
#endif /* ifndef GEDWALK_H_ */

/*
 * InfluenceMaximizationGTest.cpp
 *
 *  Created on: 17.10.2018
 *      Author: Eugenio Angriman
 */

#include <gtest/gtest.h>

#include "../../auxiliary/Log.h"
#include "../../auxiliary/PrioQueue.h"
#include "../../auxiliary/Random.h"
#include "../../centrality/Closeness.h"
#include "../../distance/CommuteTimeDistance.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../InfluenceMaximization.h"

namespace NetworKit {

class InfluenceMaximizationGTest : public testing::Test {};
TEST_F(InfluenceMaximizationGTest, testExpectedValue) {
	const count n = 20;
	const double p = 0.2;

	const count simulations = 1e6;
	const node v = 0;
	node best;
	double bestScore = 0.;

	node bestCCd;
	double bestScoreCCd = 0.;

	count seed = 100;
	// Generates a strogly connected graph
	for (count seed = 100; seed < seed + 1; seed += 2) {
		Aux::Random::setSeed(seed, false);
		Graph G(n, true, true);
		for (node u = 0; u < n - 1; ++u) {
			G.addEdge(u, u + 1, Aux::Random::real(0, 1));
		}
		// Closing loop
		G.addEdge(n - 1, 0, Aux::Random::real(0, 1));
		for (node u = 0; u < n; ++u) {
			for (node v = 0; v < n; ++v) {
				if (v != u && Aux::Random::real(0, 1) <= p) {
					G.addEdge(u, v, Aux::Random::real(0, 1));
				}
			}
		}

		INFO("Number of edges = ", G.numberOfEdges());
		Aux::PrioQueue<double, node> simPQ(n);
		Aux::PrioQueue<double, node> compPQ(n);
		InfluenceMaximization infl(G);
		infl.computeDiameter();
		for (node u = 0; u < n; ++u) {
			auto result = infl.performSimulation(u, simulations);
			double curResult = std::accumulate(result.begin(), result.end(), 0.) - 1.;
			simPQ.insert(-curResult, u);
			auto upperBound = infl.computeInflProb(u);
			double computedResult =
			    std::accumulate(upperBound.begin(), upperBound.end(), 0.);
			compPQ.insert(-computedResult, u);
		}

		for (node u = 0; u < n; ++u) {
			auto simElem = simPQ.peekMin(u);
			simElem.first *= -1.;
			auto compElem = compPQ.peekMin(u);
			compElem.first *= -1.;
			INFO("Simulated: ", simElem, " --- Upper Bound: ", compElem);
		}
	}
}

} /* namespace NetworKit */

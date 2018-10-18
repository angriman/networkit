/*
 * InfluenceMaximizationGTest.cpp
 *
 *  Created on: 17.10.2018
 *      Author: Eugenio Angriman
 */

#include <gtest/gtest.h>

#include "../../auxiliary/Log.h"
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

	const count simulations = 3e6;
	const node v = 0;
	node best;
	double bestScore = 0.;

	node bestCCd;
	double bestScoreCCd = 0.;

	count seed = 100;
	for (count seed = 100; seed < 200; seed += 2) {
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
		INFO("Total edge weight = ", G.totalEdgeWeight());
		INFO("Top Closeness ", cc.ranking()[0]);
		CommuteTimeDistance ccd(G, 0.001);
		ccd.run();
		//		InfluenceMaximization infl(G);
		//		std::vector<double> inflScores;
		std::vector<double> cfccScores;
		for (node u = 0; u < n; ++u) {
			//			auto result = infl.performSimulation(u, simulations);
			//			double curResult = std::accumulate(result.begin(), result.end(),
			// 0.); 			inflScores.push_back(curResult); 			if (curResult >
			// bestScore) { 				bestScore = curResult; 				best = u;
			//			}
			double curCCd = 0.;
			for (node v = 0; v < n; ++v) {
				curCCd += ccd.distance(u, v);
			}
			curCCd = 1. / curCCd;
			cfccScores.push_back(curCCd);
			if (curCCd > bestScoreCCd) {
				bestScoreCCd = curCCd;
				bestCCd = u;
			}
			//	INFO("Done simulation with ", u + 1, " nodes");
		}
		//		EXPECT_EQ(best, bestCCd);
		//		INFO("Highest value: ", best, " with score ", bestScore);
		INFO("Best node ccd: ", bestCCd, " with score ", bestScoreCCd);
		//		std::sort(inflScores.begin(), inflScores.end(),
		// std::greater<double>()); 		std::sort(cfccScores.begin(),
		// cfccScores.end(), std::greater<double>()); 		INFO(inflScores);
		// INFO(cfccScores); 	INFO("Computed Values:\n", infl.computeInflProb(v));
		//	INFO("Simulation result:\n", infl.performSimulation(v, simulations));
	}
}

} /* namespace NetworKit */

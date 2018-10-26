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
#include "../../centrality/WeightedDegreeCentrality.h"
#include "../../distance/CommuteTimeDistance.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../io/EdgeListReader.h"
#include "../IndependentCascade.h"
#include "../InfluenceMaximization.h"
#include "../LinearThreshold.h"
#include "../TopCentrality.h"

namespace NetworKit {

class InfluenceMaximizationGTest : public testing::Test {};

TEST_F(InfluenceMaximizationGTest, testTopCentrality) {
	Graph G(5, true, true);
	G.addEdge(0, 1, .5);
	G.addEdge(0, 2, .5);
	G.addEdge(0, 3, .5);
	G.addEdge(0, 4, .5);
	G.addEdge(2, 3, .5);
	G.addEdge(2, 1, .5);
	G.addEdge(4, 3, .5);
	G.addEdge(4, 1, .5);

	TopCentrality topc(G, 3);
	topc.run();
	topc.getInfluencers();
}

TEST_F(InfluenceMaximizationGTest, testAlgorithm) {
	//	const count n = 10;
	//	const double p = .2;
	//	Graph G(n, true, true);
	//	for (count i = 0; i < n; ++i) {
	//		for (count j = 0; j < n; ++j) {
	//			if (i != j && Aux::Random::probability() <= p) {
	//				G.addEdge(i, j, Aux::Random::probability());
	//			}
	//		}
	//	}

	Aux::Random::setSeed(1, false);
	EdgeListReader reader(' ', 0, "#", false, true);
	Graph G1 = reader.read("/home/angriman/graphs/facebook_combined");
	Graph G(G1.numberOfNodes(), true, true);
	G1.forEdges(
	    [&](node u, node v) { G.addEdge(u, v, Aux::Random::probability()); });

	InfluenceMaximization infl(G);
	infl.run();
}

TEST_F(InfluenceMaximizationGTest, testExpectedValueSimpleGraph) {
	Graph G(4, true, true);
	G.addEdge(0, 1, .5);
	G.addEdge(0, 2, .5);
	G.addEdge(1, 2, .5);
	G.addEdge(2, 1, .5);
	G.addEdge(1, 3, .5);
	G.addEdge(2, 3, .5);

	InfluenceMaximization infl(G);
	infl.computeDiameter();

	IndependentCascade IC(G);
	IC.performSimulation({0}, 1e7);
	auto simResult = IC.getAverage();
	auto expResult = infl.computeInflProb(0);

	INFO("Simulation: ", simResult);
	INFO("Computed:   ", expResult);
}

TEST_F(InfluenceMaximizationGTest, testExpectedValueBiggerGraph) {
	const count n = 20;
	const double p = .2;
	Graph G(n, true, true);
	for (count i = 0; i < n; ++i) {
		for (count j = 0; j < n; ++j) {
			if (i != j && Aux::Random::probability() <= p) {
				G.addEdge(i, j, Aux::Random::probability());
			}
		}
	}

	InfluenceMaximization infl(G);
	IndependentCascade IC(G);
	for (node v = 0; v < n; ++v) {
		IC.performSimulation({v}, 1e7);
		INFO("Simulation:   ", IC.getAverage());
		INFO("Computed val: ", infl.computeInflProb(v));
	}
}

TEST_F(InfluenceMaximizationGTest, testTree) {
	Graph g(10, true, true);

	g.addEdge(0, 1, .5);
	g.addEdge(0, 2, .5);
	g.addEdge(0, 3, .5);
	g.addEdge(1, 4, .5);
	g.addEdge(1, 5, .5);
	g.addEdge(2, 6, .5);
	g.addEdge(2, 7, .5);
	g.addEdge(2, 8, .5);
	g.addEdge(3, 9, .5);

	InfluenceMaximization infl(g);
	IndependentCascade IC(g);
	IC.performSimulation({0}, 1e8);
	INFO("Simulation: ", IC.getAverage());
	INFO("Computed:   ", infl.computeInflProb(0));
}

TEST_F(InfluenceMaximizationGTest, testDag) {
	Graph g(10, true, true);

	g.addEdge(0, 1, .5);
	g.addEdge(0, 2, .5);
	g.addEdge(0, 3, .5);
	g.addEdge(1, 4, .5);
	g.addEdge(1, 5, .5);
	g.addEdge(1, 6, .5);
	g.addEdge(2, 6, .5);
	g.addEdge(2, 7, .5);
	g.addEdge(2, 8, .5);
	g.addEdge(3, 9, .5);
	g.addEdge(3, 8, .5);

	InfluenceMaximization infl(g);
	IndependentCascade IC(g);
	IC.performSimulation({0}, 1e8);
	INFO("Simulation: ", IC.getAverage());
	INFO("Computed:   ", infl.computeInflProb(0));
}

TEST_F(InfluenceMaximizationGTest, testExpectedValue) {
	const count n = 90;
	const double p = 0.01;

	const count simulations = 5e4;
	const node v = 0;
	node best;
	double bestScore = 0.;

	node bestCCd;
	double bestScoreCCd = 0.;

	count seed = 100;
	// Generates a strogly connected graph
	for (count seed = 100; seed <= 100; seed += 2) {
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
		IndependentCascade IC(G);
		infl.computeDiameter();
		for (node u = 0; u < n; ++u) {
			IC.performSimulation({u}, simulations);
			auto curResult = IC.getAverage();
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
			INFO("Simulated: ", simElem, " --- Approximation: ", compElem);
		}

		IC.performSimulation({simPQ.peekMin(0).second}, simulations);
		INFO(IC.getAverage());
		INFO(infl.computeInflProb(simPQ.peekMin(0).second));
	}
}

TEST_F(InfluenceMaximizationGTest, testWeightedDegreIn) {
	Graph G(4, true, true);
	G.addEdge(1, 0, 0.4);
	G.addEdge(2, 0, 0.2);
	G.addEdge(3, 0, 0.1);
	G.addEdge(0, 1, 0.2);
	G.addEdge(1, 3, 0.2);

	EXPECT_NEAR(G.weightedDegreeIn(0), 0.7, 1e-6);
	EXPECT_NEAR(G.weightedDegreeIn(1), 0.2, 1e-6);
	EXPECT_NEAR(G.weightedDegreeIn(2), 0.0, 1e-6);
	EXPECT_NEAR(G.weightedDegreeIn(3), 0.2, 1e-6);
}

TEST_F(InfluenceMaximizationGTest, testWeightedDegreCentrality) {
	Graph G(4, true, true);
	G.addEdge(1, 0, 0.4);
	G.addEdge(2, 0, 0.2);
	G.addEdge(3, 0, 0.1);
	G.addEdge(0, 1, 0.2);
	G.addEdge(1, 3, 0.2);

	std::vector<double> scores(4, .1); // = {.5, .9, .2, .4};
	WeightedDegreeCentrality wdc(G, scores);
	wdc.run();
	auto ranking = wdc.ranking();
	for (count i = 0; i < ranking.size() - 1; ++i) {
		EXPECT_TRUE(ranking[i].second >= ranking[i + 1].second);
	}
	for (auto r : ranking) {
		double weightedDeg = 0.;
		G.forNeighborsOf(r.first, [&](node u) {
			weightedDeg += G.weight(r.first, u) * scores[u];
		});
		EXPECT_EQ(weightedDeg, r.second);
	}
}

} /* namespace NetworKit */

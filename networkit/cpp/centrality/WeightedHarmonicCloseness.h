/*
 * WeightedHarmonicCloseness.h
 *
 * Created on: 24.02.2018
 * 		 Author: Eugenio Angriman
 */

#ifndef WEIGHTEDHARMONICCLOSENESS_H_
#define WEIGHTEDHARMONICCLOSENESS_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class WeightedHarmonicCloseness : public NetworKit::Centrality {
public:
	/**
	 * Constructs the WeghtedHarmonicCloseness class for the given Graph @a G. If
	 * the closeness scores should be normalized, then set @a normalized to
	 * <code>true</code>. The run() method takes O(nm) time, where n is the number
	 * of nodes and m is the number of edges of the graph.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>false</code> if scores should
	 * not be normalized into an interval of [0, 1]. Normalization only for
	 * unweighted graphs.
	 *
	 */
	WeightedHarmonicCloseness(const Graph &G,
	                          const std::vector<double> &nodeWeights,
	                          const std::vector<bool> &inGroup,
	                          bool normalized = true);

	/**
	 * Computes the harmonic closeness centrality on the graph passed in
	 * constructor.
	 */
	void run() override;

	/*
	 * Returns the maximum possible harmonic closeness centrality that a node can
	 * have in a star graph with the same amount of nodes.
	 */
	double maximum() override;

private:
	const std::vector<double> &nodeWeights;
	const std::vector<bool> &inGroup;
	const count n;
};
} // namespace NetworKit

#endif

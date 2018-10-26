/* GroupCloseness.cpp
 *
 *	Created on: 18.6.2018
 *			Author: Eugenio Angriman
 */

#include <cmath>
#include <omp.h>

#include "GedWalk.h"

namespace NetworKit {

GedWalk::GedWalk(Graph &G, const count k, const double alpha)
    : G(G), n(G.upperNodeIdBound()), k(k),
      alpha(alpha <= 0. ? 1. / (1. + (double)maxDegree()) : alpha),
      alphas(initAlphas()), q(n),
      gMat(CSRMatrix::adjacencyMatrix(G).transpose()) {
	if (k == 0 || k >= n) {
		throw std::runtime_error("Error: k should be between 1 and "
		                         "n-1.");
	}
	gMat.sort();
	G.sortEdges();
}

void GedWalk::init() {
	hasRun = false;
	oldScore = newScore = 0.f;

	inGroup.assign(n, false);
	isExact.assign(n, false);
	group.reserve(k);
	score.assign(n, 0.f);

	pathsIn.resize(maxLevel, std::vector<count>(n, 0));
	pathsOut.resize(maxLevel, std::vector<count>(n, 0));

	// Each node has a 0-path
	if (G.isDirected()) {
		std::fill(pathsIn[0].begin(), pathsIn[0].end(), 1);
	}
	std::fill(pathsOut[0].begin(), pathsOut[0].end(), 1);
}

void GedWalk::computeUpperBound() {
	count contrib;
	for (uint8_t i = 1; i < maxLevel; ++i) {
#pragma omp parallel for private(contrib)
		for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
			contrib = 0;
			gMat.forNonZeroElementsInRow(u, [&](node v, edgeweight eid) {
				if (G.isDirected()) {
					contrib += pathsIn[i - 1][v];
				}
			});

			if (G.isDirected()) {
				pathsIn[i][u] += contrib;
			}

			contrib = 0;
			G.forNeighborsOf(u, [&](node v) { contrib += pathsOut[i - 1][v]; });
			pathsOut[i][u] += contrib;
		}
	}

	std::vector<std::vector<count>> maxPaths(maxLevel, std::vector<count>(n, 0));
	// Path length
	for (uint8_t l = 1; l < maxLevel; ++l) {
		count nPaths;
		uint8_t l1;
#pragma omp for private(nPaths, l1)
		for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
			nPaths = 0;
			// l1 is the suffix. Paths ending in u will have a suffix length = 0.
			for (l1 = 0; l1 <= l; ++l1) {
				if (G.isDirected()) {
					nPaths += pathsIn[l - l1][u] * pathsOut[l1][u];
				} else {
					nPaths += pathsOut[l - l1][u] * pathsOut[l1][u];
				}
			}
			maxPaths[l][u] = nPaths;
			score[u] += alphas[l - 1] * nPaths;
		}
	}
	G.forNodes([&](node u) { q.insert(-score[u], u); });
}

void GedWalk::computePaths() {
	count pIn, pOut;
	for (uint8_t i = 1; i < maxLevel; ++i) {
#pragma omp parallel for private(pIn, pOut)
		for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
			pIn = pOut = 0;
			gMat.forNonZeroElementsInRow(u, [&](node v, edgeweight eid) {
				if (inGroup[v]) {
					pIn += pathsOut[i - 1][v];
				} else {
					pOut += pathsOut[i - 1][v];
				}
				pIn += pathsIn[i - 1][v];
			});
			pathsIn[i][u] = pIn;
			pathsOut[i][u] = pOut;
		}
	}
}

double GedWalk::computeScore() {
	double result = 0.;
#pragma omp parallel for reduction(+ : result)
	for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
		double contrib = 0.;
		for (uint8_t i = 1; i < maxLevel; ++i) {
			contrib += (double)pathsIn[i][u] * alphas[i - 1];
		}
		result += contrib;
	}
	return result;
}

void GedWalk::computeMarginalGain(const node &z) {
	// For simplicity
	inGroup[z] = true;
	pathsIn[0][z] = 1;
	pathsOut[0][z] = 0;
	computePaths();
	// Restoring coherence
	inGroup[z] = false;
	pathsIn[0][z] = 0;
	pathsOut[0][z] = 1;
	double newScore = computeScore();
	score[z] = newScore - oldScore;
}

void GedWalk::run() {
	init();

	computeUpperBound();
	std::fill(pathsIn[0].begin(), pathsIn[0].end(), 0);

	node z;
	double margGain;
	std::pair<double, node> top;

	// TODO parallelize this loop
	while (group.size() < k) {
		top = q.extractMin();
		margGain = -top.first;
		z = top.second;
		if (!isExact[z]) {
			computeMarginalGain(z);
			q.insert(-score[z], z);
			isExact[z] = true;
			continue;
		}

		newScore = margGain + oldScore;
		oldScore = newScore;

		inGroup[z] = true;
		group.push_back(z);
		pathsIn[0][z] = 1;
		pathsOut[0][z] = 0;
		std::fill(isExact.begin(), isExact.end(), false);
	}

	hasRun = true;
}

} // namespace NetworKit

/*
 * SpanningEdgeCentrality.cpp
 *
 *  Created on: 29.07.2015
 *      Author: henningm
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/SpanningEdgeCentrality.hpp>

#include <fstream>
#include <sstream>

#include <omp.h>

namespace NetworKit {

SpanningEdgeCentrality::SpanningEdgeCentrality(const Graph &G, double tol)
    : Centrality(G, false, true), tol(tol), lamg(1e-5) {
    // prepare LAMG
    CSRMatrix matrix = CSRMatrix::laplacianMatrix(G);
    Aux::Timer t;
    t.start();
    lamg.setupConnected(matrix);
    t.stop();

    setupTime = t.elapsedMilliseconds();

    DEBUG("done setting up SpanningEdgeCentrality");
}

void SpanningEdgeCentrality::run() {
    count n = G.numberOfNodes();
    scoreData.clear();
    scoreData.resize(G.upperEdgeIdBound(), 0.0);

    // set up solution vector and status
    Vector solution(n);

    Vector rhs(n, 0.0);
    Vector zeroVector(n, 0.0);

    // solve for each edge
    G.forEdges([&](node u, node v, edgeweight weight, edgeid e) {
        // set up right-hand side
        rhs[u] = +1.0;
        rhs[v] = -1.0;
        TRACE("before solve for ", u, " and ", v);

        solution = zeroVector;

        lamg.solve(rhs, solution);
        double diff = solution[u] - solution[v];
        scoreData[e] = weight * std::fabs(diff);
        rhs[u] = 0.0;
        rhs[v] = 0.0;
    });

    hasRun = true;
}

uint64_t SpanningEdgeCentrality::getSetupTime() const {
    return setupTime;
}

void SpanningEdgeCentrality::runApproximation() {
    const count n = G.numberOfNodes();
    double epsilon2 = tol * tol;
    const count k = std::ceil(std::log2(n)) / epsilon2;
    double randTab[3] = {1 / std::sqrt(k), -1 / std::sqrt(k)};
    Vector solution(n);
    scoreData.clear();
    scoreData.resize(G.upperEdgeIdBound(), 0.0);

    for (index i = 0; i < k; ++i) {
        Vector rhs(n, 0.0);

        // rhs(v) = \sum_e=1 ^m q(e) * sqrt(w(e)) * B(e, v)
        //        = +/- q(e) * sqrt(w(e))
        G.forEdges([&](node u, node v, edgeweight weight) {
            double r = randTab[Aux::Random::integer(1)] * std::sqrt(weight);

            if (u < v) {
                rhs[u] += r;
                rhs[v] -= r;
            } else {
                rhs[u] -= r;
                rhs[v] += r;
            }
        });

        lamg.solve(rhs, solution);

        G.forEdges([&](node u, node v, edgeweight weight, edgeid e) {
            double diff = solution[u] - solution[v];
            scoreData[e] += weight * diff * diff;
        });
    }

    hasRun = true;
}

void SpanningEdgeCentrality::runParallelApproximation() {
    const count n = G.numberOfNodes();
    double epsilon2 = tol * tol;
    const count k = std::ceil(std::log2(n)) / epsilon2;
    double randTab[3] = {1 / std::sqrt(k), -1 / std::sqrt(k)};
    std::vector<Vector> solutions(k, Vector(n));
    std::vector<Vector> rhs(k, Vector(n));
    scoreData.clear();
    scoreData.resize(G.upperEdgeIdBound(), 0.0);

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(k); ++i) {
        // rhs(v) = \sum_e=1 ^m q(e) * sqrt(w(e)) * B(e, v)
        //        = +/- q(e) * sqrt(w(e))
        G.forEdges([&](node u, node v, edgeweight weight) {
            double r = randTab[Aux::Random::integer(1)] * std::sqrt(weight);

            if (u < v) {
                rhs[i][u] += r;
                rhs[i][v] -= r;
            } else {
                rhs[i][u] -= r;
                rhs[i][v] += r;
            }
        });
    }

    lamg.parallelSolve(rhs, solutions);

    for (index i = 0; i < k; ++i) {
        G.parallelForEdges([&](node u, node v, edgeweight weight, edgeid e) {
            double diff = solutions[i][u] - solutions[i][v];
            scoreData[e] += weight * diff * diff;
        });
    }

    hasRun = true;
}

uint64_t SpanningEdgeCentrality::runApproximationAndWriteVectors(std::string_view) {
    WARN("SpanningEdgeCentrality::runApproximationAndWriteVectors should not be used and will be "
         "deprecated in the future.");
    Aux::Timer t;
    t.start();
    runApproximation();
    t.stop();
    return t.elapsedMilliseconds();
}

double SpanningEdgeCentrality::runForEdge(node u, node v) {
    count n = G.numberOfNodes();

    // set up solution vector and status
    Vector solution(n, 0.0);
    Vector rhs(n, 0.0);

    // set up right-hand side
    rhs[u] = +1.0;
    rhs[v] = -1.0;
    TRACE("before solve for ", u, " and ", v);

    lamg.solve(rhs, solution);
    return G.weight(u, v) * std::fabs(solution[u] - solution[v]);
}

} // namespace NetworKit

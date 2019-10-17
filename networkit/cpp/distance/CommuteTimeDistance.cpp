/*
 * CommuteTimeDistance.cpp
 *
 *  Created on: 29.07.2015
 *      Author: henning
 */

#include <math.h>
#include <omp.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/distance/CommuteTimeDistance.hpp>

namespace NetworKit {

CommuteTimeDistance::CommuteTimeDistance(const Graph& G, double tol): Algorithm(), G(&G), tol(tol), lamg(1e-5) {
    // main purpose of method: preparing LAMG

    // construct matrix from graph
    CSRMatrix matrix = CSRMatrix::laplacianMatrix(G);

    // run LAMG setup and measure the time
    Aux::Timer t;
    t.start();
    lamg.setupConnected(matrix);
    t.stop();
    setupTime = t.elapsedMilliseconds();

    DEBUG("done setting up commute time distance");
    hasRun = false;
}

void CommuteTimeDistance::run() {
    count n = G->numberOfNodes();
    distances.resize(n);
    G->forNodes([&](node v){
        distances[v].resize(n, 0.0);
    });

    // temp vector for resetting the solution state
    Vector zeroVector(n, 0.0);

    // set up solution vector and right-hand side
    Vector solution = zeroVector;
    Vector rhs = zeroVector;

    // solve for each pair of nodes
    G->forNodePairs([&](node u, node v){
        // set up right-hand side
        rhs[u] = +1.0;
        rhs[v] = -1.0;
        TRACE("before solve for ", u, " and ", v);

        solution = zeroVector;

        lamg.solve(rhs, solution);
        double diff = fabs(solution[u] - solution[v]);
        distances[u][v] = diff;
        distances[v][u] = diff;
        rhs[u] = 0.0;
        rhs[v] = 0.0;
    });
    exactly = true;
    hasRun = true;
}

uint64_t CommuteTimeDistance::getSetupTime() const {
    return setupTime;
}

void CommuteTimeDistance::runApproximation() {
    count n = G->numberOfNodes();

    // init approximation parameters
    double epsilon2 = tol * tol;
    k = ceil(log2(n)) / epsilon2;

    // entries of random projection matrix
    double randTab[2] = {1 / sqrt(k), -1 / sqrt(k)};

    solutions.clear();
    solutions.resize(k, Vector(n));

    for (index i = 0; i < k; ++i) {
        Vector rhs(n, 0.0);

        // matrix vector product of q
        // rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
        //        = +/- q(e)
        G->forEdges([&](node u, node v) {
            double r = randTab[Aux::Random::integer(1)];

            if (u < v) {
                rhs[u] += r;
                rhs[v] -= r;
            } else {
                rhs[u] -= r;
                rhs[v] += r;
            }
        });

        lamg.solve(rhs, solutions[i]);
    }
    exactly = false;
    hasRun = true;
}

void CommuteTimeDistance::runParallelApproximation() {
    count n = G->numberOfNodes();

    // init approximation parameters
    double epsilon2 = tol * tol;
    k = ceil(log2(n)) / epsilon2;

    // entries of random projection matrix
    double randTab[3] = {1 / sqrt(k), -1 / sqrt(k)};

    solutions.clear();
    solutions.resize(k, Vector(n));
    std::vector<Vector> rhs(k, Vector(n));

    INFO("Number k of iterations: ", k);
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(k); ++i) {
        G->forEdges([&](node u, node v) {
            double r = randTab[Aux::Random::integer(1)];

            if (u < v) {
                rhs[i][u] += r;
                rhs[i][v] -= r;
            } else {
                rhs[i][u] -= r;
                rhs[i][v] += r;
            }
        });
    }
    INFO("Starting the solve phase");
    lamg.parallelSolve(rhs, solutions);
    INFO("Done with the solve phase");

    exactly = false;
    hasRun = true;
}

double CommuteTimeDistance::distance(node u, node v) {
    assureFinished();

    // compute volume
    double volG = 0.0;
    if (! G->isWeighted()) {
        volG = 2.0 * G->numberOfEdges();
    }
    else {
        volG = 2.0 * G->totalEdgeWeight();
    }

    if (exactly) {
        return sqrt(distances[u][v] * volG);
    } else {
        double dist = 0.0;
        for (index i = 0; i < k; ++i) {
            double diff = solutions[i][u] - solutions[i][v];
            dist += diff * diff;
        }
        return sqrt(dist * volG);
    }
}

double CommuteTimeDistance::runSinglePair(node u, node v) {
    count n = G->numberOfNodes();
    double dist = 0.0;

    // set up solution vector and status
    Vector solution(n);

    Vector rhs(n, 0.0);
    Vector zeroVector(n, 0.0);
    rhs[u] = +1.0;
    rhs[v] = -1.0;
    // set up right-hand side
    solution = zeroVector;
    lamg.solve(rhs, solution);
    double diff = solution[u] - solution[v];
    dist = fabs(diff);
    return sqrt(dist* G->numberOfEdges());
}

double CommuteTimeDistance::effectiveResistanceSinglePair(node u, node v) {
    count n = G.numberOfNodes();

    // set up solution vector and status
    Vector solution(n);

    Vector rhs(n, 0.0);
    Vector zeroVector(n, 0.0);
    rhs[u] = +1.0;
    rhs[v] = -1.0;
    // set up right-hand side
    solution = zeroVector;
    lamg.solve(rhs, solution);
    double diff = solution[u] - solution[v];
    return std::abs(diff);
}

std::vector<double> CommuteTimeDistance::effectiveResistanceSingleSource(node u) {
    const count n = G.numberOfNodes();
    std::vector<double> result(n);
    Vector solution(n);
    Vector rhs(n, 0.0);
    Vector zeroVector(n, 0.0);

    Aux::Timer timer;
    timer.start();
    G.forNodes([&](const node v) {
        rhs = zeroVector;
        solution = zeroVector;
        rhs[u] = +1.0;
        rhs[v] = -1.0;
        lamg.solve(rhs, solution);
        result[v] = std::abs(solution[u] - solution[v]);
    });

    timer.stop();
    elapsedMilliseconds = timer.elapsedMilliseconds();
    return result;
}

std::vector<double> CommuteTimeDistance::effectiveResistanceSingleSourceParallel(node u) {
    count n = G.numberOfNodes();
    std::vector<double> result(n);

    // Solution vectors: one per thread
    std::vector<Vector> solutions(omp_get_max_threads(), Vector(n));

    // Right hand side vectors: one per thread
    std::vector<Vector> rhss(omp_get_max_threads(), Vector(n, 0.0));

    // Zero vector, used to reset;
    Vector zeroVector(n, 0.0);

    Aux::Timer timer;
    timer.start();
    for (count i = 0; i <= n / static_cast<count>(omp_get_max_threads()); ++i) {
        // Index of the next vertex to process
        const index base = i * omp_get_max_threads();

#pragma omp parallel
        {
            // Each thread solves a linear system from `base` to `base + #threads - 1`
            const index v = base + omp_get_thread_num();
            if (v < n) {
                // Reset solution and rhs vector of the current thread
                solutions[omp_get_thread_num()] = zeroVector;
                rhss[omp_get_thread_num()] = zeroVector;

                // Set up system to compute the effective resistance from `u` to `idx`
                rhss[omp_get_thread_num()][u] = +1.0;
                rhss[omp_get_thread_num()][v] = -1.0;
            }
        }

        if (base + omp_get_max_threads() < n) {
            // All threads can be employed
            lamg.parallelSolve(rhss, solutions);
        } else {
            // Last iteration: some threads cannot be employed.
            // Resize rhss and solutions to the number of vertices left to be processed.
            rhss.resize(n - base + 1);
            solutions.resize(rhss.size());
            lamg.parallelSolve(rhss, solutions);
        }

        // Store the results
        for (omp_index idx = 0; idx < omp_get_max_threads(); ++idx)
            if (base + idx < n) {
                result[base + idx] = std::abs(solutions[idx][u] - solutions[idx][base + idx]);
            }
    }
    timer.stop();
    elapsedMilliseconds = timer.elapsedMilliseconds();
    return result;
>>>>>>> Create ApproxEffectiveResistance class and run test
}

double CommuteTimeDistance::runSingleSource(node u) {
    count n = G->numberOfNodes();
    double dist = 0.0;
    double sum = 0.0;
    Vector zeroVector(n, 0.0);
    // set up solution vector and status
    std::vector<Vector> rhs(n, Vector(n));
    std::vector<Vector> solution(n, Vector(n));
    G->forNodes([&](node i){
        rhs[i] = zeroVector;
        solution[i] = zeroVector;
        rhs[i][u] = +1.0;
        if (i != u) {
            rhs[i][i] = -1.0;
        } else {
            rhs[i][0] = -1.0;
        }
    });

    INFO("rhs.size() = ", rhs.size());
    INFO("solutions.size() = ", solution.size());
    lamg.parallelSolve(rhs, solution);
    G->forNodes([&](node i){
        if (i != u) {
            double diff = solution[i][u] - solution[i][i];
            dist = fabs(diff);
            sum += sqrt(dist);
        }
    });
    return sum * sqrt(G->numberOfEdges());
}

} // namespace NetworKit

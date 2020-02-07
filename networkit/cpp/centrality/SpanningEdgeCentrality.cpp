/*
 * SpanningEdgeCentrality.cpp
 *
 *  Created on: 29.07.2015
 *      Author: henningm
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/SpanningEdgeCentrality.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>

#include <fstream>
#include <random>
#include <sstream>

#include <omp.h>

namespace NetworKit {

SpanningEdgeCentrality::SpanningEdgeCentrality(const Graph &G, double tol)
    : Centrality(G), tol(tol), lamg(1e-5) {
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
    scoreData.resize(G.numberOfEdges(), 0.0);

    // set up solution vector and status
    Vector solution(n);

    Vector rhs(n, 0.0);
    Vector zeroVector(n, 0.0);

    // solve for each edge
    G.forEdges([&](node u, node v, edgeid e) {
        // set up right-hand side
        rhs[u] = +1.0;
        rhs[v] = -1.0;
        TRACE("before solve for ", u, " and ", v);

        solution = zeroVector;

        lamg.solve(rhs, solution);
        double diff = solution[u] - solution[v];
        scoreData[e] = fabs(diff); // TODO: check unweighted, fix weighted case, fix edge IDs!
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
    const count m = G.numberOfEdges();
    double epsilon2 = tol * tol;
    const count k = ceil(log2(n)) / epsilon2;
    double randTab[3] = {1 / sqrt(k), -1 / sqrt(k)};
    Vector solution(n);
    scoreData.clear();
    scoreData.resize(m, 0.0);

    for (index i = 0; i < k; ++i) {
        Vector rhs(n, 0.0);

        // rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
        //        = +/- q(e)
        G.forEdges([&](node u, node v) {
            double r = randTab[Aux::Random::integer(1)];

            if (u < v) {
                rhs[u] += r;
                rhs[v] -= r;
            } else {
                rhs[u] -= r;
                rhs[v] += r;
            }
        });

        lamg.solve(rhs, solution);

        G.forEdges([&](node u, node v, edgeid e) {
            double diff = solution[u] - solution[v];
            scoreData[e] += diff * diff; // TODO: fix weighted case!
        });
    }

    hasRun = true;
}

  
  
// Approximation algo to compute diagonal of pseudoinverse explicitly.
// Based on random vectors ( Hutchinson estimator ).
// Input: # of samples, tolerance.
  std::vector<double> SpanningEdgeCentrality::computeDiagonalHadaEst(int k, double tol) {
    if ((k % 4)!= 0)
      throw std::runtime_error("Error. \n Number of samples must be a multiple of 4.");
    
    const auto n = G.numberOfNodes();
    int numbits = ceil(log(n)/log(2));
    const int next_pow = pow(2, numbits);
    numbits++;
    // create the hadamard binary representation
    std::vector<std::vector<double>> Hadabin(next_pow, std::vector<double>(numbits));
    for (int l = 0; l < next_pow; l++) {
      int i = 0;
      int n = l;
      while (n > 0) { 
	Hadabin[l][i] = n % 2;
	n = n / 2;
	i++;
      }
    }
        
    const omp_index n_threads = omp_get_max_threads();
    std::vector<Vector> rhss(n_threads, Vector(G.upperNodeIdBound()));
    std::vector<Vector> solutions(n_threads, Vector(G.upperNodeIdBound()));
    Lamg<CSRMatrix> solver(tol);
    solver.setupConnected(CSRMatrix::laplacianMatrix(G));

    std::vector<std::vector<double>> diags(n_threads, std::vector<double>(G.upperNodeIdBound(), 0));

    for (omp_index i = 0; i < static_cast<omp_index>(k); i += n_threads) {
#pragma omp parallel
      {
	// create Hadamard vectors from formula: H(k,j) = (-1)^( dot(Hadabin(j,:), Hadabin(k,:)));
	auto &rhs = rhss[omp_get_thread_num()];
	Vector t  = Vector(Hadabin[omp_get_thread_num()]);
	for (uint64_t i = 0; i < G.upperNodeIdBound(); i++) {
	  rhs[i] = pow(-1, Vector::innerProduct(t,Vector(Hadabin[i]))); 
	}       
      }
      solver.parallelSolve(rhss, solutions);
      
#pragma omp parallel
      {
	auto &diag = diags[omp_get_thread_num()];
	auto &solution = solutions[omp_get_thread_num()];
	G.forNodes([&](const node u) { diag[u] += solution[u] * solution[u]; });
      }
    }
    
    G.parallelForNodes([&diags, k](const node u) {
        for (size_t i = 1; i < diags.size(); ++i) {
	  diags[0][u] += diags[i][u];
        }
        diags[0][u] /= static_cast<double>(k);
      });
    
    return diags[0];
}



  


  
  
// Approximation algo to compute diagonal of pseudoinverse explicitly.
// Based on random vectors ( Hutchinson estimator ).
// Input: # of samples, tolerance.
  std::vector<double> SpanningEdgeCentrality::computeDiagonalRandomEst(int k, double tol) {
    const auto n = G.numberOfNodes();

    count n_half = std::floor(n/2);
    
    const omp_index n_threads = omp_get_max_threads();
    
    std::vector<Vector> rhss(n_threads);
    std::vector<Vector> solutions(n_threads, Vector(G.upperNodeIdBound()));
    Lamg<CSRMatrix> solver(tol);
    solver.setupConnected(CSRMatrix::laplacianMatrix(G));

    std::vector<std::vector<double>> diags(n_threads, std::vector<double>(G.upperNodeIdBound(), 0));

    for (omp_index i = 0; i < static_cast<omp_index>(k); i += n_threads) {
#pragma omp parallel
      {

	auto &rhs = rhss[omp_get_thread_num()];
	std::vector<double> in(G.upperNodeIdBound(), 1.0);
	std::fill(in.begin() + n_half, in.end(), -1.0);
	auto rng = std::default_random_engine {}; 
	std::shuffle(in.begin(), in.end(), rng);

	rhs = Vector(in);

	if ((n % 2) == 1) {
	double rsum = 0.0;
	rhs.forElements([&](const double &value) {
	    rsum += value;
	  });
	rsum = rsum/n;
	rhs.forElements([&](double &value) { value -= rsum; });
	}
      }
      solver.parallelSolve(rhss, solutions);
      
#pragma omp parallel
        {
            auto &diag = diags[omp_get_thread_num()];
            auto &solution = solutions[omp_get_thread_num()];
            G.forNodes([&](const node u) { diag[u] += solution[u] * solution[u]; });
        }
    }

    G.parallelForNodes([&diags, k](const node u) {
        for (size_t i = 1; i < diags.size(); ++i) {
            diags[0][u] += diags[i][u];
        }
        diags[0][u] /= static_cast<double>(k);
    });

    return diags[0];
}


  
std::vector<double> SpanningEdgeCentrality::computeDiagonal(double epsilon, double tol) {
    const auto n = G.numberOfNodes();
    const count k = std::ceil(log2(n)) / (epsilon * epsilon);
    static constexpr double p = 0.5;
    const double proj = 1. / std::sqrt(static_cast<double>(k));
    const omp_index n_threads = omp_get_max_threads();

    std::vector<std::mt19937_64> generators(n_threads, Aux::Random::getURNG());
    std::vector<std::uniform_real_distribution<double>> distributions(
        n_threads, std::uniform_real_distribution<double>(0, 1));
    std::vector<Vector> rhss(n_threads, Vector(G.upperNodeIdBound()));
    std::vector<Vector> solutions(n_threads, Vector(G.upperNodeIdBound()));
    Lamg<CSRMatrix> solver(tol);
    solver.setupConnected(CSRMatrix::laplacianMatrix(G));

    std::vector<std::vector<double>> diags(n_threads, std::vector<double>(G.upperNodeIdBound(), 0));

    for (omp_index i = 0; i < static_cast<omp_index>(k); i += n_threads) {
#pragma omp parallel
        {
            auto &gen = generators[omp_get_thread_num()];
            auto &distr = distributions[omp_get_thread_num()];
            auto &rhs = rhss[omp_get_thread_num()];

            rhs.forElements([](double &elem) { elem = 0; });
            G.forEdges([&](const node u, const node v) {
                const auto r = distr(gen) >= p ? proj : -proj;
                if (u < v) {
                    rhs[u] += r;
                    rhs[v] -= r;
                } else {
                    rhs[u] -= r;
                    rhs[v] += r;
                }
            });
        }
        solver.parallelSolve(rhss, solutions);

#pragma omp parallel
        {
            auto &diag = diags[omp_get_thread_num()];
            auto &solution = solutions[omp_get_thread_num()];
            G.forNodes([&](const node u) { diag[u] += solution[u] * solution[u]; });
        }
    }

    G.parallelForNodes([&diags, k](const node u) {
        for (size_t i = 1; i < diags.size(); ++i) {
            diags[0][u] += diags[i][u];
        }
        diags[0][u] /= static_cast<double>(k);
    });

    return diags[0];
}

void SpanningEdgeCentrality::runParallelApproximation() {
    const count n = G.numberOfNodes();
    const count m = G.numberOfEdges();
    double epsilon2 = tol * tol;
    const count k = ceil(log2(n)) / epsilon2;
    double randTab[3] = {1 / sqrt(k), -1 / sqrt(k)};
    std::vector<Vector> solutions(k, Vector(n));
    std::vector<Vector> rhs(k, Vector(n));
    scoreData.clear();
    scoreData.resize(m, 0.0);

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(k); ++i) {
        // rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
        //        = +/- q(e)
        G.forEdges([&](node u, node v) {
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

    lamg.parallelSolve(rhs, solutions);

    for (index i = 0; i < k; ++i) {
        G.parallelForEdges([&](node u, node v, edgeid e) {
            double diff = solutions[i][u] - solutions[i][v];
            scoreData[e] += diff * diff; // TODO: fix weighted case!
        });
    }

    hasRun = true;
}

uint64_t SpanningEdgeCentrality::runApproximationAndWriteVectors(const std::string &) {
    Aux::Timer t;
    const count n = G.numberOfNodes();
    const count m = G.numberOfEdges();
    const double epsilon2 = tol * tol;
    const count k = ceil(log(n)) / epsilon2;
    double randTab[3] = {1 / sqrt(k), -1 / sqrt(k)};
    Vector solution(n);
    scoreData.clear();
    scoreData.resize(m, 0.0);

    t.start();
    for (index i = 0; i < k; ++i) {
        Vector rhs(n, 0.0);

        // rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
        //        = +/- q(e)
        G.forEdges([&](node u, node v) {
            double r = randTab[Aux::Random::integer(1)];
            if (u < v) {
                rhs[u] += r;
                rhs[v] -= r;
            } else {
                rhs[u] -= r;
                rhs[v] += r;
            }
        });

        lamg.solve(rhs, solution);

        G.forEdges([&](node u, node v, edgeid e) {
            double diff = solution[u] - solution[v];
            scoreData[e] += diff * diff; // TODO: fix weighted case!
        });
    }
    t.stop();
    hasRun = true;

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
    return fabs(solution[u] - solution[v]); // TODO: fix weighted case, fix edge IDs!
}

} // namespace NetworKit

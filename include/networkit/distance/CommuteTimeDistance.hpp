/*
 * CommuteTimeDistance.hpp
 *
 *  Created on: 12.04.2016
 *      Author: ebergamini
 */

#ifndef NETWORKIT_DISTANCE_COMMUTE_TIME_DISTANCE_HPP_
#define NETWORKIT_DISTANCE_COMMUTE_TIME_DISTANCE_HPP_

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 *
 * CommuteTimeDistance edge centrality.
 *
 */
class CommuteTimeDistance final : public Algorithm {

public:
    /**
     * Constructs the CommuteTimeDistance class for the given Graph @a G.
     * @param G The graph.
     * @param tol The tolerance used for the approximation
     */
    CommuteTimeDistance(const Graph& G, double tol = 0.1, double lamgTol = 1e-5);

    /**
     * Destructor.
     */
    virtual ~CommuteTimeDistance() = default;


    /**
     * Computes ECTD exactly.
     */
    void run() override;

    /**
     * Computes approximation by projection.
     */
    void runApproximation();

    /**
     * Computes approximation by projection, in parallel.
     */
    void runParallelApproximation();

    /**
     * @return The elapsed time to setup the solver in milliseconds.
     */
    uint64_t getSetupTime() const;

    /**
     * Returns the commute time distance between node @a u and node @a v.
     * @return commute time distance between the two nodes. Needs to call run() or runApproximation() first.
     */
    double distance(node u, node v);

    /**
     * Returns the commute time distance between node @a u and node @a v.
     * This method does not need the initial preprocessing (no need to call the run() method).
     * @return commute time distance between the two nodes.
     */

    double runSinglePair(node u, node v);

    /**
     * Returns the sum of the distances from node @a u.
     * This method does not need the initial preprocessing.
     * @return commute sum of the distances from the node.
     */
    double runSingleSource(node u);

    double effectiveResistanceSinglePair(node u, node v);

    std::vector<double> effectiveResistanceSingleSourceParallel(node u);
    std::vector<double> effectiveResistanceSingleSource(node u);

    std::vector<double> getDiagonal(node root) {
        const auto r = effectiveResistanceSingleSourceParallel(root);
        Aux::Timer timer;
        timer.start();
        std::vector<double> diagonal(G.upperNodeIdBound());
        Vector rhs(G.upperNodeIdBound()), result(G.upperNodeIdBound());
        rhs[root] = 1.0;
        G.parallelForNodes(
            [&](const node u) { rhs[u] -= 1.0 / static_cast<double>(G.upperNodeIdBound()); });
        lamg.solve(rhs, result);
        G.parallelForNodes(
            [&](const node u) { diagonal[u] = r[u] - result[root] + 2 * result[u]; });
        timer.stop();
        elapsedMilliseconds += timer.elapsedMilliseconds();
        return diagonal;
    }

    count getElapsedMilliseconds() const {
        return elapsedMilliseconds;
    }

protected:
    const Graph* G;
    double tol, lamgTol;
    Lamg<CSRMatrix> lamg;
    uint64_t setupTime;
    std::vector<std::vector<double>> distances;
    std::vector<Vector> solutions;
    bool exactly;
    count k;
    count elapsedMilliseconds;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_COMMUTE_TIME_DISTANCE_HPP_

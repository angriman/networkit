/*
 * ConjugateGradient.h
 *
 *  Created on: 15.06.2014
 *      Author: Daniel Hoske and Michael Wegner
 */

#ifndef NETWORKIT_NUMERICS_CONJUGATE_GRADIENT_HPP_
#define NETWORKIT_NUMERICS_CONJUGATE_GRADIENT_HPP_

#include <cstdint>
#include <utility>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/numerics/LinearSolver.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 * Implementation of Conjugate Gradient.
 */
template <class Matrix, class Preconditioner>
class ConjugateGradient : public LinearSolver<Matrix> {
public:
    ConjugateGradient(double tolerance = 1e-5)
        : LinearSolver<Matrix>(tolerance), matrix(Matrix()) {}

    void setup(const Matrix &matrix) {
        this->matrix = matrix;
        precond = Preconditioner(matrix);
    }

    void setupConnected(const Matrix &matrix) {
        this->matrix = matrix;
        precond = Preconditioner(matrix);
    }

    /**
     * Solves the linear system \f$Ax = b\f$ using the conjugate gradient method
     * with a given preconditioner and with initial value \f$(0, \dots, 0)^T\f$.
     * We the return the solution \f$x\f$. The solution \f$x\f$ fulfils
     * \f$\frac{\Vert Ax - b\Vert}{\Vert b \Vert} \leq relative\_residual\f$ if the
     * algorithm has converged.
     *
     * Obviously, @a A needs to have the same number of rows as @a b and
     * @a status.residual must be nonnegative. You may also request that the algorithm
     * does not run for more than @a status.max_iters iterations.
     */
    SolverStatus solve(const Vector &rhs, Vector &result, count maxConvergenceTime = 5 * 60 * 1000,
                       count maxIterations = std::numeric_limits<count>::max());

    /**
     * Solves the linear systems in parallel.
     * @param rhs
     * @param results
     * @param maxConvergenceTime
     * @param maxIterations
     */
    void parallelSolve(const std::vector<Vector> &rhs, std::vector<Vector> &results,
                       count maxConvergenceTime = 5 * 60 * 1000,
                       count maxIterations = std::numeric_limits<count>::max());

private:
    Matrix matrix;
    Preconditioner precond;
};

template <class Matrix, class Preconditioner>
SolverStatus ConjugateGradient<Matrix, Preconditioner>::solve(const Vector &rhs, Vector &result,
                                                              count, count maxIterations) {
    assert(matrix.numberOfRows() == rhs.getDimension());

    // Absolute residual to achieve
    const double sqr_desired_residual =
        this->tolerance * this->tolerance * (rhs.length() * rhs.length());

    // Main loop. See:
    // http://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithm
    Vector residual_dir(result.getDimension(), 0); // = rhs - matrix*result;
    // TODO exchange rowIdx <-> i
#pragma omp parallel for
    for (omp_index rowIdx = 0; rowIdx < static_cast<omp_index>(matrix.numberOfRows()); ++rowIdx) {
        for (index i = matrix.rowIdxAt(rowIdx); i < matrix.rowIdxAt(rowIdx + 1); ++i) {
            residual_dir[rowIdx] += matrix.nonZerosAt(i) * result[matrix.columnIdxAt(i)];
        }
        residual_dir[rowIdx] = -residual_dir[rowIdx] + rhs[rowIdx];
    }

    Vector conjugate_dir = precond.rhs(residual_dir);
    double sqr_residual = Vector::innerProduct(residual_dir, residual_dir);
    double sqr_residual_precond = Vector::innerProduct(residual_dir, conjugate_dir);

    count niters = 0;

    const auto n = result.getDimension();
    Vector tmp(n, 0), residual_precond(n, 0);

    while (sqr_residual > sqr_desired_residual) {
        if (niters == maxIterations) {
            break;
        }
        ++niters;

        //        tmp = matrix * conjugate_dir;
        //        double step = 0; sqr_residual_precond / Vector::innerProduct(conjugate_dir, tmp);

        // Computing denominator of \alpha_k (i.e., p_k^T A p_k)
        double denominator = 0;
#pragma omp parallel for reduction(+ : denominator)
        for (omp_index row = 0; row < static_cast<omp_index>(n); ++row) {
            // Multiply A(row:) with conjugate_dir (i.e, get the resulting vector of A * p_k)
            double vecEntry = 0;
            matrix.forNonZeroElementsInRow(
                row, [&vecEntry, &conjugate_dir](const index colIdx, const double elem) {
                    vecEntry += elem * conjugate_dir[colIdx];
                });
            tmp[row] = vecEntry;
            denominator += conjugate_dir[row] * vecEntry;
        }

        // \alpha_k
        double step = sqr_residual_precond / denominator;

        //        result += step * conjugate_dir;
        //        residual_dir -= step * tmp;
        //        sqr_residual = Vector::innerProduct(residual_dir, residual_dir);

        // Update x_k, r_k and the overall residual
        sqr_residual = 0;
#pragma omp parallel for reduction(+ : sqr_residual)
        for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
            result[i] += step * conjugate_dir[i]; // x_{k+1} <-- x_k + \alpha_k * p_k
            residual_dir[i] -= step * tmp[i];     // r_{k + 1} <-- r_k - \alpha_k * A * p_k
            sqr_residual += residual_dir[i] * residual_dir[i];
        }

        if (sqr_residual <= sqr_desired_residual)
            break;

        residual_precond = precond.rhs(residual_dir);

        const double new_sqr_residual_precond = Vector::innerProduct(residual_dir, residual_precond);
        //        conjugate_dir =
        //            (new_sqr_residual_precond / sqr_residual_precond) * conjugate_dir +
        //            residual_precond;

        const double residual_fract = new_sqr_residual_precond / sqr_residual_precond;
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
            conjugate_dir[i] = residual_fract * conjugate_dir[i] + residual_precond[i];
        }

        sqr_residual_precond = new_sqr_residual_precond;
    }

    SolverStatus status;
    status.numIters = niters;
    status.residual = (rhs - matrix * result).length();
    status.converged = status.residual / rhs.length() <= this->tolerance;

    return status;
}

template <class Matrix, class Preconditioner>
void ConjugateGradient<Matrix, Preconditioner>::parallelSolve(const std::vector<Vector> &rhs,
                                                              std::vector<Vector> &results,
                                                              count maxConvergenceTime,
                                                              count maxIterations) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(rhs.size()); ++i) {
        this->solve(rhs[i], results[i], maxConvergenceTime, maxIterations);
    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_CONJUGATE_GRADIENT_HPP_

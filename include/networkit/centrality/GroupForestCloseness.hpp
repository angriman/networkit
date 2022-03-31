#ifndef NETWORKIT_CENTRALITY_GROUP_FOREST_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_GROUP_FOREST_CLOSENESS_HPP_

#include <networkit/algebraic/CSRGeneralMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/numerics/ConjugateGradient.hpp>
#include <networkit/numerics/Preconditioner/DiagonalPreconditioner.hpp>

#include <vector>

namespace NetworKit {
class GroupForestCloseness final : public Algorithm {
public:
    GroupForestCloseness(const Graph &G, node root, count k, double epsilon = 0.1);

    /**
     * Run the algorithm.
     */
    void run() override;

private:
    const Graph *G;       // Input graph
    const node root;      // Root node of the augmented graph
    const count k;        // Size of the group
    const double epsilon; // Approximation

    std::vector<node> group, nodesOutsideGroup;
    std::vector<count> neighborsInGroup;

    // Right-hand side and solution vectors used by the CG solver to compute the diagonal of the
    // pseudoinverse of graph H.
    Vector rhs, sol;

    ConjugateGradient<CSRGeneralMatrix<double>, DiagonalPreconditioner> cg;

    // Computes the first node to go in the group (i.e., the one with highest forest centrality).
    node getTopNode() const;
    void addNodeToGroup(index i);
    double computeFarnessWithNode(node u);
    double computeTraceOfPseudoInv(const Graph &H);
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_FOREST_CLOSENESS_HPP_

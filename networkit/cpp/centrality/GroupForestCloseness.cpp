#include <algorithm>
#include <cassert>
#include <limits>

#include <networkit/centrality/ApproxElectricalCloseness.hpp>
#include <networkit/centrality/ForestCentrality.hpp>
#include <networkit/centrality/GroupForestCloseness.hpp>

namespace NetworKit {
GroupForestCloseness::GroupForestCloseness(const Graph &G, node root, count k, double epsilon)
    : G(&G), root(root), k(k), epsilon(epsilon) {
    if (G.isWeighted())
        WARN("GroupForestCloseness ignores edge weights.");

    if (G.isDirected())
        throw std::runtime_error("The input graph must be undirected.");

    if (G.numberOfNodes() != G.upperNodeIdBound())
        throw std::runtime_error("The graph must be compact. Either make sure that the graph has "
                                 "no removed nodes, or use GraphTools::getCompactedGraph.");

    if (G.degree(root) != G.numberOfNodes() - 1)
        throw std::runtime_error("The input graph is not an augmented graph. Create an augmented "
                                 "graph with GraphTools::createAugmentedGraph.");
}

node GroupForestCloseness::getTopNode() const {
    ForestCentrality fc(*G, root, epsilon);
    fc.run();
    const auto &diag = fc.getDiagonal();

    // Pick non-root node with mininum forest farness
    node topNode = none;
    double minVal = std::numeric_limits<edgeweight>::max();
    G->forNodes([&](node u) {
        if (u != root && diag[u] < minVal) {
            minVal = diag[u];
            topNode = u;
        }
    });

    return topNode;
}

void GroupForestCloseness::addNodeToGroup(index i) {
    assert(i < nodesOutsideGroup.size());
    group.push_back(nodesOutsideGroup[i]);
    std::swap(nodesOutsideGroup[i], nodesOutsideGroup.back());
    nodesOutsideGroup.pop_back();
}

double GroupForestCloseness::computeFarnessWithNode(node u) {
    // Create new graph H as described in the paper.
    Graph H(*G);
    const auto n = G->numberOfNodes();
    std::fill(neighborsInGroup.begin(), neighborsInGroup.end(), 0);

    for (node v : group) {
        G->forNeighborsOf(v, [&](node neighbor) { ++neighborsInGroup[neighbor]; });
        H.removeNode(v);
    }
}

double GroupForestCloseness::computeTraceOfPseudoInv(const Graph &H) {
    ApproxElectricalCloseness apxElClos(H, root, epsilon);
    apxElClos.run();
    const auto &diagH = apxElClos.getDiagonal();
    const auto n = H.numberOfNodes();

    rhs.setDimension(n);
    sol.setDimension(n);
    rhs[root] = 1.;
    H.parallelForNodes([&](node u) { rhs[u] -= 1.0 / static_cast<double>(n); });

    cg.setTolerance(epsilon / (10.0 * std::sqrt((double)(n * H.numberOfEdges()) * std::log(n))));
    const auto L = CSRGeneralMatrix<double>::laplacianMatrix(H);
    cg.setupConnected(L);
    cg.solve(rhs, sol);

    const auto trace = H.parallelSumForNodes([&](node u) {
            return static_cast<double>
            });
}

void GroupForestCloseness::run() {
    group.clear();
    group.reserve(k);

    nodesOutsideGroup.assign(G->nodeRange().begin(), std::prev(G->nodeRange().end(), 1));
    neighborsInGroup.assign(G->numberOfNodes(), 0);

    // Add node with highest forest closeness to the group
    addNodeToGroup(getTopNode());

    while (group.size() < k) {
        index bestCandidateIdx = nodesOutsideGroup.size();
        double lowestFarness = std::numeric_limits<double>::max();

        for (index i = 0; i < nodesOutsideGroup.size(); ++i) {
            const node u = nodesOutsideGroup[i];
            const auto farnessWithU = computeFarnessWithNode(u);
            if (farnessWithU < lowestFarness) {
                lowestFarness = farnessWithU;
                bestCandidateIdx = i;
            }
        }

        assert(bestCandidateIdx < nodesOutsideGroup.size());
        addNodeToGroup(bestCandidateIdx);
    }

    hasRun = true;
}

} // namespace NetworKit

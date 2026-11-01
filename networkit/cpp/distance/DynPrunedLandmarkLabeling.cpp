#include <algorithm>
#include <cassert>
#include <cmath>
#include <queue>
#include <stdexcept>
#include <vector>

#include <networkit/distance/DynPrunedLandmarkLabeling.hpp>
#include <networkit/distance/PrunedLandmarkLabeling.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

DynPrunedLandmarkLabeling::DynPrunedLandmarkLabeling(const Graph &G)
    : PrunedLandmarkLabeling(G), ranks(G.upperNodeIdBound(), none) {
    for (node rank = 0; rank < nodesSortedByDegreeDesc.size(); ++rank) {
        ranks[nodesSortedByDegreeDesc[rank]] = rank;
    }
}

void DynPrunedLandmarkLabeling::update(GraphEvent e) {
    if (e.type == GraphEvent::EDGE_ADDITION) {
        addEdge(e.u, e.v);
    } else if (e.type == GraphEvent::EDGE_REMOVAL) {
        removeEdge(e.u, e.v);
    } else {
        throw std::runtime_error("Unsupported graph event " + e.toString());
    }
}

void DynPrunedLandmarkLabeling::sortUpdatedLabels(bool reverse) {
    for (const node u : updatedNodes) {
        auto &labelsU = reverse ? labelsIn[u] : labelsOut[u];
        if (labelsU.size() < 2) // Nothing to be sorted
            continue;

        // Insertion sort step: insert the new label in the right position within the vector.
        auto lb =
            std::lower_bound(labelsU.begin(), labelsU.end() - 2, labelsU.back(),
                             [](const auto &l1, const auto &l2) { return l1.node_ < l2.node_; });

        if (lb->node_ == labelsU.back().node_ && lb->distance_ > labelsU.back().distance_) {
            // Overwrite label
            *lb = labelsU.back();
            labelsU.pop_back();
        } else if (lb->node_ > labelsU.back().node_) {
            // Insert new label
            index i = labelsU.size() - 1, start = (lb - labelsU.begin());
            const auto newLabel = labelsU.back();
            while (i > start) {
                labelsU[i] = labelsU[i - 1];
                --i;
            }
            labelsU[i] = newLabel;
        }
    }
}

void DynPrunedLandmarkLabeling::prunedBFS(node k, node startNode, count bfsLevel, bool reverse) {
    const node root = nodesSortedByDegreeDesc[k];

    updatedNodes.clear();
    std::fill(visited.begin(), visited.end(), false);
    visited[startNode] = true;

    std::queue<node> q0, q1;
    q0.push(startNode);

    auto visitNeighbor = [&visited = visited, &q1](node v) -> void {
        if (visited[v])
            return;
        visited[v] = true;
        q1.push(v);
    };

    do {
        do {
            const node u = q0.front();
            q0.pop();

            if (reverse) {
                if (queryImpl(u, root, k) <= bfsLevel)
                    continue;
            } else {
                if (queryImpl(root, u, k) <= bfsLevel)
                    continue;
            }

            updatedNodes.push_back(u);

            if (reverse) {
                labelsIn[u].emplace_back(k, bfsLevel);
                G->forInNeighborsOf(u, visitNeighbor);
            } else {
                labelsOut[u].emplace_back(k, bfsLevel);
                G->forNeighborsOf(u, visitNeighbor);
            }
        } while (!q0.empty());

        ++bfsLevel;
        std::swap(q0, q1);
    } while (!q0.empty());

    sortUpdatedLabels(reverse);
}

void DynPrunedLandmarkLabeling::addEdge(node u, node v) {
    const bool directed = G->isDirected();
    const auto &labelsU = labelsOut[u];
    const auto &labelsV = directed ? labelsIn[v] : labelsOut[v];

    labelsUCopy.resize(labelsU.size());
    labelsVCopy.resize(labelsV.size());
    std::copy(labelsU.begin(), labelsU.end(), labelsUCopy.begin());
    std::copy(labelsV.begin(), labelsV.end(), labelsVCopy.begin());

    auto iterLabelsU = labelsUCopy.begin(), iterLabelsV = labelsVCopy.begin();
    const auto iterLabelsUEnd = labelsUCopy.end(), iterLabelsVEnd = labelsVCopy.end();

    // The searches must follow the order of the hubs. Otherwise, a search can reach a node with a
    // smaller rank than its root before the pair is covered, and the root becomes a hub of that
    // node. Removing edges relies on hubs always having a smaller rank than the nodes they cover.
    while (iterLabelsU != iterLabelsUEnd && iterLabelsV != iterLabelsVEnd) {
        if (iterLabelsU->node_ < iterLabelsV->node_) {
            prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, false);
            ++iterLabelsU;
        } else if (iterLabelsU->node_ > iterLabelsV->node_) {
            prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, directed);
            ++iterLabelsV;
        } else {
            if (directed) {
                prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, false);
                prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, true);
            } else if (iterLabelsU->distance_ + 1 < iterLabelsV->distance_) {
                prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, false);
            } else {
                prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, false);
            }

            ++iterLabelsU;
            ++iterLabelsV;
        }
    }

    while (iterLabelsU != iterLabelsUEnd) {
        prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, false);
        ++iterLabelsU;
    }

    while (iterLabelsV != iterLabelsVEnd) {
        prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, directed);
        ++iterLabelsV;
    }
}

std::vector<count> DynPrunedLandmarkLabeling::computeDistances(node source, bool reverse, node u,
                                                               node v) const {
    std::vector<count> dist(G->upperNodeIdBound(), infDist);
    std::queue<node> q;
    dist[source] = 0;
    q.push(source);

    while (!q.empty()) {
        const node w = q.front();
        q.pop();

        const auto relax = [&](node z) -> void {
            if (dist[z] == infDist) {
                dist[z] = dist[w] + 1;
                q.push(z);
            }
        };

        if (reverse) {
            G->forInNeighborsOf(w, relax);
            if (w == v) {
                relax(u);
            }
        } else {
            G->forNeighborsOf(w, relax);
            if (w == u) {
                relax(v);
            }
            if (!G->isDirected() && w == v) {
                relax(u);
            }
        }
    }

    return dist;
}

node DynPrunedLandmarkLabeling::smallestHub(node u, node v, count dist) const {
    auto [iterLabelsU, iterLabelsUEnd] = getSourceLabelsIterators(u, G->isDirected());
    auto [iterLabelsV, iterLabelsVEnd] = getSourceLabelsIterators(v);

    while (iterLabelsU != iterLabelsUEnd && iterLabelsV != iterLabelsVEnd) {
        if (iterLabelsU->node_ < iterLabelsV->node_) {
            ++iterLabelsU;
        } else if (iterLabelsV->node_ < iterLabelsU->node_) {
            ++iterLabelsV;
        } else {
            if (iterLabelsU->distance_ + iterLabelsV->distance_ == dist) {
                return iterLabelsU->node_;
            }
            ++iterLabelsU;
            ++iterLabelsV;
        }
    }

    return none;
}

std::vector<node> DynPrunedLandmarkLabeling::detectAffectedNodes(
    node source, node target, bool reverse, const std::vector<count> &distToTarget,
    const std::vector<count> &newDistToTarget, const std::vector<count> &distToSource,
    std::vector<bool> &isAffected) {
    std::vector<node> affected{source};
    std::fill(visited.begin(), visited.end(), false);
    visited[source] = true;
    isAffected[source] = true;

    const auto visitNeighbor = [&](node u) -> void {
        if (visited[u]) {
            return;
        }

        bool isAffectedU = newDistToTarget[u] != distToTarget[u];
        if (!isAffectedU) {
            const node rank = reverse ? smallestHub(u, target, distToTarget[u])
                                      : smallestHub(target, u, distToTarget[u]);
            assert(rank != none);
            const node hub = nodesSortedByDegreeDesc[rank];
            isAffectedU =
                isAffected[hub]
                || ((hub == u || hub == target) && distToTarget[u] == distToSource[u] + 1);
        }

        if (isAffectedU) {
            visited[u] = true;
            isAffected[u] = true;
            affected.push_back(u);
        }
    };

    // The affected vector doubles as the queue of the BFS.
    for (index i = 0; i < affected.size(); ++i) {
        if (reverse) {
            G->forInNeighborsOf(affected[i], visitNeighbor);
        } else {
            G->forNeighborsOf(affected[i], visitNeighbor);
        }
    }

    return affected;
}

void DynPrunedLandmarkLabeling::removeHubs(std::vector<std::vector<Label>> &labels,
                                           const std::vector<node> &nodes,
                                           const std::vector<bool> &isHub) {
    for (const node u : nodes) {
        auto &labelsU = labels[u];
        labelsU.erase(std::remove_if(labelsU.begin(), labelsU.end(),
                                     [&](const Label &label) -> bool {
                                         const node hub = nodesSortedByDegreeDesc[label.node_];
                                         return hub != u && isHub[hub];
                                     }),
                      labelsU.end());
    }
}

void DynPrunedLandmarkLabeling::insertLabel(std::vector<Label> &labelsU, node hub, count dist) {
    const auto position =
        std::lower_bound(labelsU.begin(), labelsU.end(), hub,
                         [](const Label &label, node h) -> bool { return label.node_ < h; });
    labelsU.emplace(position, hub, dist);
}

void DynPrunedLandmarkLabeling::restoreCover(node root, const std::vector<bool> &isTarget,
                                             bool reverse, bool ordered) {
    // For directed graphs, labelsIn stores the distances to the hubs and labelsOut the distances
    // from the hubs.
    auto &labelsToHub = G->isDirected() ? labelsIn : labelsOut;

    std::fill(visited.begin(), visited.end(), false);
    visited[root] = true;

    std::queue<node> q0, q1;
    q0.push(root);
    count level = 0;

    const auto visitNeighbor = [&](node w) -> void {
        if (!visited[w]) {
            visited[w] = true;
            q1.push(w);
        }
    };

    do {
        do {
            const node v = q0.front();
            q0.pop();

            if (ordered && ranks[v] < ranks[root]) {
                continue;
            }

            if (isTarget[v] && level < (reverse ? queryImpl(v, root) : queryImpl(root, v))) {
                if (ranks[v] < ranks[root]) {
                    insertLabel(reverse ? labelsOut[root] : labelsToHub[root], ranks[v], level);
                } else {
                    insertLabel(reverse ? labelsToHub[v] : labelsOut[v], ranks[root], level);
                }
            }

            if (reverse) {
                G->forInNeighborsOf(v, visitNeighbor);
            } else {
                G->forNeighborsOf(v, visitNeighbor);
            }
        } while (!q0.empty());

        ++level;
        std::swap(q0, q1);
    } while (!q0.empty());
}

void DynPrunedLandmarkLabeling::removeEdge(node x, node y) {
    // Self-loops and parallel edges do not affect distances.
    if (x == y || G->hasEdge(x, y)) {
        return;
    }

    const bool directed = G->isDirected();

    const auto toX = computeDistances(x, directed, x, y);
    const auto toY = computeDistances(y, directed, x, y);
    const auto newToY = computeDistances(y, directed);
    const auto fromX = directed ? computeDistances(x, false, x, y) : toX;
    const auto fromY = directed ? computeDistances(y, false, x, y) : toY;
    const auto newFromX = computeDistances(x, false);

    // Phase 1: nodes with a shortest path through the removed edge that is induced by the labels.
    std::vector<bool> isAffectedX(G->upperNodeIdBound(), false);
    std::vector<bool> isAffectedY(G->upperNodeIdBound(), false);
    const auto affectedX = detectAffectedNodes(x, y, directed, toY, newToY, toX, isAffectedX);
    const auto affectedY = detectAffectedNodes(y, x, false, fromX, newFromX, fromY, isAffectedY);

    // Phase 2: only labels between nodes of the two sets can be outdated.
    removeHubs(directed ? labelsIn : labelsOut, affectedX, isAffectedY);
    removeHubs(labelsOut, affectedY, isAffectedX);

    // If the endpoints are disconnected, so are all pairs of nodes of the two sets. This does not
    // hold for directed graphs, where the two sets can be connected by paths that avoid x and y.
    if (!directed && newToY[x] == infDist) {
        return;
    }

    // Phase 3: Order-Restore is faster unless both sets are large (roughly n / log n according to
    // the experiments in the paper).
    const count n = G->numberOfNodes();
    const count minAffected = std::min(affectedX.size(), affectedY.size());
    if (static_cast<double>(minAffected) <= static_cast<double>(n) / std::log2(n)) {
        std::vector<node> roots(affectedX);
        roots.insert(roots.end(), affectedY.begin(), affectedY.end());
        std::sort(roots.begin(), roots.end(),
                  [this](node u, node v) -> bool { return ranks[u] < ranks[v]; });
        roots.erase(std::unique(roots.begin(), roots.end()), roots.end());

        for (const node root : roots) {
            if (isAffectedX[root]) {
                restoreCover(root, isAffectedY, /*reverse=*/false, /*ordered=*/true);
            }
            if (isAffectedY[root]) {
                restoreCover(root, isAffectedX, /*reverse=*/directed, /*ordered=*/true);
            }
        }
    } else if (affectedX.size() <= affectedY.size()) {
        for (const node root : affectedX) {
            restoreCover(root, isAffectedY, /*reverse=*/false, /*ordered=*/false);
        }
    } else {
        for (const node root : affectedY) {
            restoreCover(root, isAffectedX, /*reverse=*/directed, /*ordered=*/false);
        }
    }
}

} // namespace NetworKit

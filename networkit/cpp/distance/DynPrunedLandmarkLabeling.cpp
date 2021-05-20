// networkit-format

#include <networkit/distance/DynPrunedLandmarkLabeling.hpp>

namespace NetworKit {

constexpr count DynPrunedLandmarkLabeling::infDist;

DynPrunedLandmarkLabeling::DynPrunedLandmarkLabeling(const Graph &G)
    : G(&G), nodesSortedByDegree(G.nodeRange().begin(), G.nodeRange().end()) {

    if (G.isDirected())
        Aux::Parallel::sort(nodesSortedByDegree.begin(), nodesSortedByDegree.end(),
                            [&G](const node u, const node v) {
                                const auto degU = G.degree(u), degV = G.degree(v);
                                if (degU != degV)
                                    return degU > degV;
                                return G.degreeIn(u) > G.degreeIn(v);
                            });
    else
        Aux::Parallel::sort(nodesSortedByDegree.begin(), nodesSortedByDegree.end(),
                            [&G](const node u, const node v) { return G.degree(u) > G.degree(v); });

    rankOfNode.resize(G.upperNodeIdBound());
    node curRank = 0;

    for (const node u : nodesSortedByDegree)
        rankOfNode[u] = curRank++;

    visited.resize(G.upperNodeIdBound());
    labelsOut.resize(G.upperNodeIdBound());
    updatedNodes.reserve(G.upperNodeIdBound());

    if (G.isDirected())
        labelsIn.resize(G.upperNodeIdBound());

    labelsUCopy.reserve(G.upperNodeIdBound());
    labelsVCopy.reserve(G.upperNodeIdBound());

#ifndef NDEBUG
    apsp = std::unique_ptr<APSP>(new APSP(G));
    apsp->run();
#endif
}

void DynPrunedLandmarkLabeling::updateAffected() {
    if (tsAffected++ == maxTS) {
        tsAffected = 1;
        std::fill(affected.begin(), affected.end(), 0);
    }
}

void DynPrunedLandmarkLabeling::prunedBFS(node root, node rootRank, bool reverse) {
    std::fill(visited.begin(), visited.end(), false);
    visited[root] = true;

    std::queue<node> q0, q1;
    q0.push(root);

    count level = 0;
    do {
        do {
            const node u = q0.front();
            q0.pop();
            if (reverse) {
                if (u != root && query(u, root) <= level)
                    continue;
            } else if (u != root && query(root, u) <= level)
                continue;

            auto visitNeighbor = [&visited = visited, &q1](node v) -> void {
                if (visited[v])
                    return;
                visited[v] = true;
                q1.push(v);
            };

            if (reverse) {
                labelsIn[u].emplace_back(rootRank, level);
                G->forInNeighborsOf(u, visitNeighbor);
                assert(apsp->getDistance(u, root) == level);
            } else {
                labelsOut[u].emplace_back(rootRank, level);
                G->forNeighborsOf(u, visitNeighbor);
                assert(apsp->getDistance(root, u) == level);
            }

        } while (!q0.empty());

        ++level;
        std::swap(q0, q1);
    } while (!q0.empty());
}

void DynPrunedLandmarkLabeling::sortUpdatedLabels() {
    for (const node u : updatedNodes) {
        auto &labelsU = G->isDirected() ? labelsIn[u] : labelsOut[u];
        if (labelsU.size() < 2)
            continue;
        const auto lb =
            std::lower_bound(labelsU.begin(), labelsU.end() - 2, labelsU.back(),
                             [](const auto &l1, const auto &l2) { return l1.node_ < l2.node_; });

        if (lb->node_ == labelsU.back().node_ && lb->distance_ > labelsU.back().distance_) {
            // Overwrite label
            *lb = std::move(labelsU.back());
            labelsU.pop_back();
        } else if ((*lb).node_ > labelsU.back().node_) {

            // Insert new label
            count i = labelsU.size() - 1, start = (lb - labelsU.begin());
            const auto newLabel = std::move(labelsU.back());
            while (i > start) {
                labelsU[i] = std::move(labelsU[i - 1]);
                --i;
            }

            labelsU[i] = std::move(newLabel);
        }
    }
}

void DynPrunedLandmarkLabeling::resumeBFS(node k, node startNode, node otherNode, count level,
                                          bool reverse) {
    const node root = nodesSortedByDegree[k];

    updatedNodes.clear();
    std::fill(visited.begin(), visited.end(), false);
    visited[startNode] = true;
    std::queue<node> q0, q1;
    q0.push(startNode);

    do {
        do {
            const node u = q0.front();
            q0.pop();
            if (reverse) {
                if (query(u, root, true, k) <= level)
                    continue;
            } else {
                if (query(root, u, true, k) <= level)
                    continue;
            }

            auto visitNeighbor = [&visited = visited, &q1](node v) -> void {
                if (visited[v])
                    return;
                visited[v] = true;
                q1.push(v);
            };

            updatedNodes.push_back(u);
            if (reverse) {
                assert(apsp->getDistance(u, root) == level);
                labelsIn[u].emplace_back(k, level);
                G->forInNeighborsOf(u, visitNeighbor);
            } else {
                assert(apsp->getDistance(root, u) == level);
                labelsOut[u].emplace_back(k, level);
                G->forNeighborsOf(u, visitNeighbor);
            }
        } while (!q0.empty());

        ++level;
        std::swap(q0, q1);
    } while (!q0.empty());

    sortUpdatedLabels();
}

void DynPrunedLandmarkLabeling::computePrevAndNewDists(node x, node y, std::vector<count> &prevDist,
                                                       std::vector<count> &newDist) {
    std::fill(prevDist.begin(), prevDist.end(), infDist);
    std::fill(newDist.begin(), newDist.end(), infDist);

    prevDist[x] = 0;
    prevDist[y] = 1;
    std::queue<node> q;
    q.push(x);
    q.push(y);

    const auto doBFS = [&](auto &distance) -> void {
        do {
            const node u = q.front();
            q.pop();
            G->forNeighborsOf(u, [&](const node v) {
                if (distance[v] == infDist) {
                    distance[v] = distance[u] + 1;
                    q.push(v);
                }
            });
        } while (!q.empty());
    };

    doBFS(prevDist);
    q.push(x);
    newDist[x] = 0;
    doBFS(newDist);
#ifndef NDEBUG
    G->forNodes([&](const auto u) {
        assert(apsp->getDistance(x, u) == prevDist[u]);
        assert(apspNew->getDistance(x, u) == newDist[u]);
    });
#endif
}

void DynPrunedLandmarkLabeling::computeAffectedNodes(node x, node y,
                                                     std::vector<node> &affectedNodesRanks,
                                                     std::vector<count> &prevDist,
                                                     std::vector<count> &newDist,
                                                     std::vector<count> &prevDist2) {

    std::fill(visited.begin(), visited.end(), false);
    visited[x] = true;

    updateAffected();
    affectedNodesRanks.clear();

    std::queue<node> q0, q1;
    q0.push(x);

    count level = 0;
    do {
        do {
            const node v = q0.front();
            q0.pop();
            affected[v] = tsAffected;
            affectedNodesRanks.push_back(rankOfNode[v]);

            auto visitNeighbor = [&](node u) {
                if (visited[u])
                    return;
                if (prevDist[u] != newDist[u]) {
                    q1.push(u);
                    visited[u] = true;
                } else {
                    const node h = smallestHub(u, y);
                    assert(apsp->getDistance(x, u) == prevDist2[u]);
                    if (affected[h] == tsAffected
                        || ((h == u || h == y) && prevDist[u] == prevDist2[u] + 1)) {
                        q1.push(u);
                        visited[u] = true;
                    }
                }
            };
            G->forNeighborsOf(v, visitNeighbor);
        } while (!q0.empty());

        ++level;
        std::swap(q0, q1);
    } while (!q0.empty());
}

void DynPrunedLandmarkLabeling::computeNewHubs(node rootRank, uint8_t tsAffectedU,
                                               uint8_t tsAffectedV, bool rootInAffectedU) {
    const node root = nodesSortedByDegree[rootRank];
    std::fill(visited.begin(), visited.end(), false);
    visited[root] = true;
    updatedNodes.clear();

    std::queue<node> q0, q1;
    q0.push(root);

    count level = 0;
    do {
        do {
            const node v = q0.front();
            q0.pop();

            const auto insertIfFarther = [&]() {
                if (level < query(root, v)) {
                    assert(level == apsp->getDistance(root, v));
                    labelsOut[v].emplace_back(rootRank, level);
                    updatedNodes.push_back(v);
                }
            };

            if (rootInAffectedU) {
                if (affected[v] == tsAffectedV)
                    insertIfFarther();
            } else if (affected[v] == tsAffectedU)
                insertIfFarther();

            G->forNeighborsOf(v, [&](const node u) {
                if (visited[u])
                    return;
                visited[u] = true;
                if (rankOfNode[u] < rootRank)
                    return;
                q1.push(u);
            });
        } while (!q0.empty());

        ++level;
        std::swap(q0, q1);
    } while (!q0.empty());

    sortUpdatedLabels();
}

void DynPrunedLandmarkLabeling::removeAffectedHubs(uint8_t tsAffectedU, uint8_t tsAffectedV) {
    auto removeAffected = [&](node xRank, uint8_t tsToTest) -> void {
        auto &labelsX = labelsOut[nodesSortedByDegree[xRank]];
        int removedLabels = 0, lowestToRemove = labelsX.size();

        for (int i = labelsX.size() - 1; i >= 0; --i) {
            if (affected[nodesSortedByDegree[labelsX[i].node_]] == tsToTest) {
                ++removedLabels;
                lowestToRemove = i;
                labelsX[i].node_ = none;
            }
        }

        for (int i = lowestToRemove + 1; i < labelsX.size(); ++i) {
            if (labelsX[i].node_ == none)
                continue;
            labelsX[lowestToRemove++] = std::move(labelsX[i]);
        }

        labelsX.resize(labelsX.size() - removedLabels);
    };

    for (const node x : affectedURanks)
        removeAffected(x, tsAffectedV);

    for (const node x : affectedVRanks)
        removeAffected(x, tsAffectedU);

#ifndef NDEBUG
    checkLabels();
#endif
}

count DynPrunedLandmarkLabeling::query(node u, node v, bool prefixal, node upperBound) const {
    assert(G->hasNode(u) && G->hasNode(v));

    auto iterLabelsU = G->isDirected() ? labelsIn[u].begin() : labelsOut[u].begin();
    const auto iterLabelsUEnd = G->isDirected() ? labelsIn[u].end() : labelsOut[u].begin();
    auto iterLabelsV = labelsOut[v].begin();
    const auto iterLabelsVEnd = labelsOut[v].end();

    count result = infDist;
    while (iterLabelsU != iterLabelsUEnd && iterLabelsV != iterLabelsVEnd) {
        if (prefixal)
            if (std::max(iterLabelsU->node_, iterLabelsV->node_) > upperBound)
                break;
        if (iterLabelsU->node_ < iterLabelsV->node_)
            ++iterLabelsU;
        else if (iterLabelsU->node_ > iterLabelsV->node_)
            ++iterLabelsV;
        else {
            assert(iterLabelsU != iterLabelsUEnd && iterLabelsV != iterLabelsVEnd);
            const node x = iterLabelsU->node_;
            result = std::min(result, iterLabelsU->distance_ + iterLabelsV->distance_);

            ++iterLabelsU;
            ++iterLabelsV;
        }
    }

    return result;
}

node DynPrunedLandmarkLabeling::smallestHub(node u, node v) const {
    auto iterLabelsU = labelsOut[u].begin(), iterLabelsV = labelsOut[v].begin();
    node smallestHubRank = none;
    count minDist = infDist;
    do {
        if (iterLabelsU->node_ < iterLabelsV->node_)
            ++iterLabelsU;
        else if (iterLabelsU->node_ > iterLabelsV->node_)
            ++iterLabelsV;
        else {
            const count newDist = iterLabelsU->distance_ + iterLabelsV->distance_;
            if (newDist < minDist) {
                minDist = newDist;
                smallestHubRank = iterLabelsU->node_;
            }
            ++iterLabelsU;
            ++iterLabelsV; } } while (iterLabelsU != labelsOut[u].end() && iterLabelsV != labelsOut[v].end());

    assert(smallestHubRank != none);
    return nodesSortedByDegree[smallestHubRank];
}

void DynPrunedLandmarkLabeling::addEdge(node u, node v) {
    assert(!G->hasEdge(u, v));
    assert(G->hasNode(u) && G->hasNode(v));
#ifndef NDEBUG
    auto G1 = *G;
    G1.addEdge(u, v);
    apsp = std::unique_ptr<APSP>(new APSP(G1));
    apsp->run();
#endif
    const auto &labelsU = labelsOut[u], &labelsV = G->isDirected() ? labelsIn[v] : labelsOut[v];

    labelsUCopy.resize(labelsU.size());
    labelsVCopy.resize(labelsV.size());
    std::copy(labelsU.begin(), labelsU.end(), labelsUCopy.begin());
    std::copy(labelsV.begin(), labelsV.end(), labelsVCopy.begin());

    auto iterLabelsU = labelsUCopy.begin(), iterLabelsV = labelsVCopy.begin();
    const auto iterLabelsUEnd = labelsUCopy.end(), iterLabelsVEnd = labelsVCopy.end();

    assert(iterLabelsU != iterLabelsUEnd);
    assert(iterLabelsV != iterLabelsVEnd);

    if (!G->isDirected())
        do {
            if (iterLabelsU->node_ < iterLabelsV->node_) {
                resumeBFS(iterLabelsU->node_, v, u, iterLabelsU->distance_ + 1);
                ++iterLabelsU;
            } else if (iterLabelsU->node_ > iterLabelsV->node_) {
                resumeBFS(iterLabelsV->node_, u, v, iterLabelsV->distance_ + 1);
                ++iterLabelsV;
            } else {
                if (iterLabelsU->distance_ + 1 < iterLabelsV->distance_)
                    resumeBFS(iterLabelsU->node_, v, u, iterLabelsU->distance_ + 1);
                else
                    resumeBFS(iterLabelsV->node_, u, v, iterLabelsV->distance_ + 1);

                ++iterLabelsU;
                ++iterLabelsV;
            }
        } while (iterLabelsU != iterLabelsUEnd && iterLabelsV != iterLabelsVEnd);

    while (iterLabelsU != iterLabelsUEnd) {
        resumeBFS(iterLabelsU->node_, v, u, iterLabelsU->distance_ + 1);
        ++iterLabelsU;
    }

    while (iterLabelsV != iterLabelsVEnd) {
        resumeBFS(iterLabelsV->node_, u, v, iterLabelsV->distance_ + 1, G->isDirected());
        ++iterLabelsV;
    }

#ifndef NDEBUG
    checkLabels();
#endif
}

void DynPrunedLandmarkLabeling::removeEdge(node u, node v) {
    assert(!G->hasEdge(u, v));
    assert(G->hasNode(u) && G->hasNode(v));
#ifndef NDEBUG
    apspNew = std::unique_ptr<APSP>(new APSP(*G));
    apspNew->run();
#endif

    prevDistX.resize(G->upperNodeIdBound());
    prevDistY.resize(G->upperNodeIdBound());
    newDistX.resize(G->upperNodeIdBound());
    newDistY.resize(G->upperNodeIdBound());
    affected.resize(G->upperNodeIdBound(), 0);
    computePrevAndNewDists(u, v, prevDistX, newDistX);
    computePrevAndNewDists(v, u, prevDistY, newDistY);

    computeAffectedNodes(u, v, affectedURanks, prevDistY, newDistY, prevDistX);
    const auto tsAffectedU = tsAffected;
    computeAffectedNodes(v, u, affectedVRanks, prevDistX, newDistX, prevDistY);
    const auto tsAffectedV = tsAffected;

#ifndef NDEBUG
    if (!G->isDirected())
        for (node x : affectedURanks)
            assert(std::find(affectedVRanks.begin(), affectedVRanks.end(), x)
                   == affectedVRanks.end());
    apsp->run();
    checkAffectedHubs(tsAffectedU, tsAffectedV);
    G->forNodes([&](const auto x) {
        node rankX = rankOfNode[x];
        if (affected[x] == tsAffectedU)
            assert(std::find(affectedURanks.begin(), affectedURanks.end(), rankX)
                   != affectedURanks.end());
        else if (affected[x] == tsAffectedV)
            assert(std::find(affectedVRanks.begin(), affectedVRanks.end(), rankX)
                   != affectedVRanks.end());
    });
#endif

    removeAffectedHubs(tsAffectedU, tsAffectedV);
    std::sort(affectedURanks.begin(), affectedURanks.end());
    std::sort(affectedVRanks.begin(), affectedVRanks.end());
    auto iterU = affectedURanks.begin(), iterV = affectedVRanks.begin();

    do {
        if (*iterU < *iterV)
            computeNewHubs(*(iterU++), tsAffectedU, tsAffectedV, true);
        else
            computeNewHubs(*(iterV++), tsAffectedU, tsAffectedV, false);
    } while (iterU != affectedURanks.end() && iterV != affectedVRanks.end());

    while (iterU != affectedURanks.end())
        computeNewHubs(*(iterU++), tsAffectedU, tsAffectedV, true);
    while (iterV != affectedVRanks.end())
        computeNewHubs(*(iterV++), tsAffectedU, tsAffectedV, false);
#ifndef NDEBUG
    checkLabels();
#endif
}

void DynPrunedLandmarkLabeling::run() {
    node rootRank = 0;
    for (const node root : nodesSortedByDegree) {
        prunedBFS(root, rootRank);
        if (G->isDirected())
            prunedBFS(root, rootRank, true);
        ++rootRank;
    }
#ifndef NDEBUG
    checkLabels();
#endif
    hasRun = true;
}

#ifndef NDEBUG
void DynPrunedLandmarkLabeling::checkLabels() const {
    auto checkLabelsVec = [&](const auto &labelsU) -> void {
        node u = 0;
        for (const auto &l : labelsU) {
            for (int i = 1; i < l.size(); ++i)
                assert(l[i - 1].node_ < l[i].node_);
            ++u;
        }
    };
    checkLabelsVec(labelsOut);
    checkLabelsVec(labelsIn);
}

void DynPrunedLandmarkLabeling::checkAffectedHubs(uint8_t tsAffectedU, uint8_t tsAffectedV) const {
    G->forNodePairs([&](const auto u, const auto v) {
        auto condition = [&]() -> bool {
            return (affected[u] == tsAffectedU && affected[v] == tsAffectedV)
                   || (affected[u] == tsAffectedV && affected[v] == tsAffectedU);
        };
        if (query(u, v) != apsp->getDistance(u, v))
            assert(condition());
    });
}
#endif

} // namespace NetworKit

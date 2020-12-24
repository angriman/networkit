// networkit-format

#include <atomic>
#include <cmath>
#include <omp.h>
#include <queue>

#include "GroupClosenessImpl.hpp"

#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {
namespace GroupClosenessDetails {

template <class Weight>
GroupClosenessImpl<Weight>::GroupClosenessImpl(const Graph &G, count k, Weight maxDist)
    : G(&G), k(k), maxDist(maxDist) {

    if (k == 0 || k >= G.numberOfNodes())
        throw std::runtime_error("Error, k must be in [0, n - 1].");

    const count n = G.upperNodeIdBound(), threads = omp_get_max_threads();
    group.reserve(k);
    distFromGroup.resize(n, infDist);
    distGlobal.resize(threads, std::vector<Weight>(n));
    lowestDist.resize(n, infDist);
    visitedGlobal.resize(threads, std::vector<bool>(n));
    farness.resize(n, 0);
    marginalGain.resize(n);

    nearerNodesGlobal.resize(threads);
    for (count i = 0; i < threads; ++i)
        nearerNodesGlobal[i].reserve(n - 1);

    candidateNodesPQ.reserve(n);
    if (G.isWeighted()) {
        dijkstraHeaps.reserve(threads);
        for (count i = 0; i < threads; ++i) {
            dijkstraHeaps.emplace_back(distGlobal[i]);
            dijkstraHeaps.back().reserve(n);
            minEdgeWeight = infDist;
            G.forEdges([&minEdgeWeight = minEdgeWeight](node, node, Weight ew) {
                minEdgeWeight = std::min(minEdgeWeight, ew);
            });
        }
    } else
        predGlobal.resize(threads, std::vector<node>(n));
}

template <class Weight>
void GroupClosenessImpl<Weight>::run() {
    group.push_back(topClosenessNode());
#ifdef NETWORKIT_SANITY_CHECKS
    checkDistFromGroup();
#endif // NETWORKIT_SANITY_CHECKS

    while (group.size() < k) {

        candidateNodesPQ.build_heap(G->nodeRange().begin(), G->nodeRange().end());
        for (node u : group)
            candidateNodesPQ.remove(u);

        group.push_back(findNodeWithHighestMarginalGain());
#ifdef NETWORKIT_SANITY_CHECKS
        checkDistFromGroup();
#endif // NETWORKIT_SANITY_CHECKS
        marginalGain[group.back()] = 0;
    }
}

template <class Weight>
node GroupClosenessImpl<Weight>::topClosenessNode() {
    computeFarnessLowerBound();
    tlx::DAryAddressableIntHeap<node, 2, Less> nodesWithLowestFarness{Less(farness)};
    nodesWithLowestFarness.build_heap(G->nodeRange().begin(), G->nodeRange().end());
    node highestClosenessNode = none;
    double lowestFarness = std::numeric_limits<double>::max();

    std::atomic<bool> stop{false};

#pragma omp parallel
    {
        while (!stop.load(std::memory_order_relaxed)) {
            node u = none;
#pragma omp critical
            {
                if (!nodesWithLowestFarness.empty()) {
                    u = nodesWithLowestFarness.extract_top();
                    if (farness[u] >= lowestFarness) {
                        stop.store(true, std::memory_order_relaxed);
                        u = none;
                    }
                } else
                    stop.store(true, std::memory_order_relaxed);
            }

            if (u == none)
                break;

            const auto ssspResult = prunedSSSPEmptyGroup(u, lowestFarness);
            marginalGain[u] = ssspResult.farness > 0 ? 1. / static_cast<double>(ssspResult.farness) : 0;
            if (!ssspResult.pruned) {
#ifdef NETWORKIT_SANITY_CHECKS
                checkTopNode(u, distGlobal[omp_get_thread_num()], ssspResult.farness);
#endif // NETWORKIT_SANITY_CHECKS

#pragma omp critical
                if (ssspResult.farness < lowestFarness) {
                    lowestFarness = ssspResult.farness;
                    highestClosenessNode = u;
                    std::swap(distFromGroup, distGlobal[omp_get_thread_num()]);
                }
            }
        }
    }

    assert(highestClosenessNode != none);
    return highestClosenessNode;
}

template <>
GroupClosenessImpl<count>::PrunedSSSPResult
GroupClosenessImpl<count>::prunedSSSPEmptyGroup(node source, double lowestFarness) {
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);
    visited[source] = true;

    auto &dist = distGlobal[omp_get_thread_num()];
    dist[source] = 0;

    auto &pred = predGlobal[omp_get_thread_num()];
    std::queue<node> q1, q2;
    q1.push(source);

    count level = 1, visitedNodes = 1, curFarness = 0;
    const count undirected = static_cast<count>(!G->isDirected());
    count farnessLowerBound = farness[source];

    do {
        count nodesAtNextLevelUB = 0, prevAtNextLevelUB = 0;
        do {
            const node u = q1.front();
            q1.pop();

            for (const node w : G->neighborRange(u)) {
                if (!visited[w]) {
                    visited[w] = true;
                    q2.push(w);
                    ++visitedNodes;
                    curFarness += level;
                    nodesAtNextLevelUB += G->degree(w) - undirected;
                    if (!G->isDirected())
                        pred[w] = u;
                } else if (!G->isDirected() && prevAtNextLevelUB > 0 && u != pred[w]) {
                    assert(u != source);
                    --prevAtNextLevelUB;
                    ++farnessLowerBound;
                    if (farnessLowerBound >= lowestFarness)
                        return {false, farnessLowerBound};
                }
            }
        } while (!q1.empty());

        assert(visitedNodes <= G->numberOfNodes());
        nodesAtNextLevelUB = std::min(nodesAtNextLevelUB, G->numberOfNodes() - visitedNodes);
        prevAtNextLevelUB = nodesAtNextLevelUB;
        farnessLowerBound =
            curFarness                         // Farness computed so far
            + nodesAtNextLevelUB * (level + 1) // Upper bound of the vertices at next level
            + (G->numberOfNodes() - visitedNodes - nodesAtNextLevelUB) * (level + 2); // Remaining

        if (farnessLowerBound >= lowestFarness)
            return {false, farnessLowerBound};

        ++level;
        std::swap(q1, q2);
    } while (!q1.empty());

    return {true, curFarness};
}

template <>
GroupClosenessImpl<edgeweight>::PrunedSSSPResult
GroupClosenessImpl<edgeweight>::prunedSSSPEmptyGroup(node source, double lowestFarness) {
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);
    visited[source] = true;

    auto &dist = distGlobal[omp_get_thread_num()];
    std::fill(dist.begin(), dist.end(), infDist);
    dist[source] = 0;

    auto &prioQ = dijkstraHeaps[omp_get_thread_num()];
    prioQ.clear();
    prioQ.push(source);

    const auto exploreNeighbors = [&](node u) -> void {
        node w;
        edgeweight ew;
        for (const auto neighborWeight : G->weightNeighborRange(u)) {
            std::tie(w, ew) = neighborWeight;
            const edgeweight newDist = dist[u] + ew;
            if (!visited[w]) {
                visited[w] = true;
                dist[w] = newDist;
                prioQ.push(w);
            } else if (newDist < dist[w]) {
                dist[w] = newDist;
                prioQ.update(w);
            }
        }
    };

    edgeweight curFarness = 0, farnessLowerBound = farness[source];
    count visitedNodes = 1;

    do {
        const node u = prioQ.extract_top();
        curFarness += dist[u];
        farnessLowerBound =
            curFarness + static_cast<edgeweight>(G->numberOfNodes() - visitedNodes) * dist[u];

        if (farnessLowerBound >= lowestFarness)
            return {false, farnessLowerBound};

        ++visitedNodes;
        exploreNeighbors(u);
    } while (!prioQ.empty());

    return {true, curFarness};
}

template <>
node GroupClosenessImpl<count>::findNodeWithHighestMarginalGain() {
    return none;
}

template <>
node GroupClosenessImpl<edgeweight>::findNodeWithHighestMarginalGain() {
    return none;
}

template <>
void GroupClosenessImpl<count>::computeFarnessLowerBound() {
    G->parallelForNodes([&farness = farness, &G = G](node u) {
        farness[u] = G->degree(u) + (G->numberOfNodes() - G->degree(u)) * 2;
    });
}

template <>
void GroupClosenessImpl<edgeweight>::computeFarnessLowerBound() {
    G->parallelForNodes([&farness = farness, &G = G, &minEdgeWeight = minEdgeWeight](node u) {
        edgeweight cheapestOutEdge = infDist;
        G->forNeighborsOf(u, [&cheapestOutEdge](node, edgeweight ew) {
            cheapestOutEdge = std::min(cheapestOutEdge, ew);
        });

        farness[u] =
            cheapestOutEdge
            + (minEdgeWeight + cheapestOutEdge) * static_cast<edgeweight>(G->numberOfNodes() - 1);
    });
}

#ifdef NETWORKIT_SANITY_CHECKS
template <>
void GroupClosenessImpl<count>::checkTopNode(node u, const std::vector<count> &dist,
                                             count computedFarness) const {
    count farnessU = 0;
    Traversal::BFSfrom(*G, u, [&dist = dist](node v, count distV) {
        assert(dist[v] == distV);
        farnessU += distV;
    });
    assert(farnessU == computedFarness);
}

template <>
void GroupClosenessImpl<edgweight>::checkTopNode(node u, const std::vector<edgweight> &dist,
                                                 edgeweight computedFarness) const {
    edgeweight farnessU = 0;
    Traversal::DijkstraFrom(*G, u, [&dist = dist](node v, edgweight distV) {
        assert(dist[v] == distV);
        farnessU += distV;
    });
    assert(farnessU == computedFarness);
}

template <>
void GroupClosenessImpl<count>::checkDistFromGroup() const {
    Traversal::BFSfrom(
        *G, group.begin(), group.end(),
        [&distFromGroup = distFromGroup](node u, count dist) { assert(distFromGroup[u] == dist); });
}

template <>
void GroupClosenessImpl<edgeweight>::checkDistFromGroup() const {
    Traversal::DijkstraFrom(*G, group.begin(), group.end(),
                            [&distFromGroup = distFromGroup](node u, edgeweight dist) {
                                assert(distFromGroup[u] == dist);
                            });
}
#endif // NETWORKIT_SANITY_CHECKS

template class GroupClosenessImpl<count>;
template class GroupClosenessImpl<edgeweight>;

} // namespace GroupClosenessDetails
} // namespace NetworKit

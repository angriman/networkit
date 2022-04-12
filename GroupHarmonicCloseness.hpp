#ifndef NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_HPP_

#include <atomic>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <queue>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/Dijkstra.hpp>
#include <networkit/graph/Graph.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {
template <class Weight>
class GroupHarmonicClosenessGeneral final : public Algorithm {
    static constexpr Weight infDist = std::numeric_limits<Weight>::std::max();

public:
    const Graph &G;
    const count k;
    std::vector<Weight> distFromGroup, curBestDist, nextDistFromGroup;
    std::vector<std::vector<Weight>> curDistGlobal;
    std::vector<count> reachableNodesInComponent, reachableNodesUB;
    std::vector<index> componentOfNode;
    std::vector<node> pred, group, nearerNodes;
    std::vector<std::vector<node>> curNearerNodesGlobal;
    std::vector<double> margGain;
    std::vector<std::vector<uint8_t>> visitedGlobal;
    std::vector<uint8_t> tsGlobal;
    Weight minEdgeWeight = infDist;

    template <class Type>
    struct Greater {
        Greater(const std::vector<Type> &margGain) : margGain(margGain) {}
        bool operator()(node x, node y) const noexcept { return margGain[x] > margGain[y]; }

    private:
        const std::vector<Type> &margGain;
    };

    template <class Type>
    struct Less {
        Less(const std::vector<Type> &dist) : dist(dist) {}
        bool operator()(node x, node y) const noexcept { return dist[x] < dist[y]; }

    private:
        const std::vector<Type> &dist;
    };

    tlx::d_ary_addressable_int_heap<node, 2, Greater<double>> candidateNodesPQ{
        Greater<double>(margGain)};
    std::vector<tlx::d_ary_addressable_int_heap<node, 2, Less<Weight>>> dijkstraHeaps;

    void incTS() {
        auto &ts = tsGlobal[omp_get_thread_num()];
        if (ts++ == std::numeric_limits<uint8_t>::std::max()) {
            ts = 1;
            auto &visited = visitedGlobal[omp_get_thread_num()];
            std::fill(visited.begin(), visited.end(), 0);
        }
    }

    void computeReachableNodesUndir();
    void computeReachableNodesDir();

    std::pair<bool, double> bfscut(node source, double bestScore);
    std::pair<bool, double> prunedBFS(node source, double bestScore);
    node findTopHCloseness();
    node findNodeWithHighestMargGain();

    double closenessUBUnd(node u) const noexcept;
    double closenessUBDir(node u) const noexcept;

    GroupHarmonicClosenessGeneral(const Graph &G, count k = 1) : G(G), k(k) {
        if (k == 0)
            throw std::runtime_error("k must be greater than zero");

        const count n = G.upperNodeIdBound();
        margGain.resize(n);
        distFromGroup.resize(n, infDist);
        curDistGlobal.resize(omp_get_max_threads(), std::vector<Weight>(n, infDist));
        curBestDist.resize(n);
        tsGlobal.resize(omp_get_max_threads(), 0);
        visitedGlobal.resize(omp_get_max_threads(), std::vector<uint8_t>(n, 0));
        candidateNodesPQ.reserve(n);
        group.reserve(k);
        nearerNodes.reserve(n - 1);
        curNearerNodesGlobal.resize(omp_get_max_threads());
        for (omp_index i = 0; i < omp_get_max_threads(); ++i)
            curNearerNodesGlobal[i].reserve(n - 1);

        componentOfNode.resize(n);

        if (G.isWeighted()) {
            dijkstraHeaps.reserve(omp_get_max_threads());
            for (omp_index i = 0; i < omp_get_max_threads(); ++i) {
                dijkstraHeaps.emplace_back(curDistGlobal[i]);
                dijkstraHeaps.back().reserve(n);
            }
            G.forEdges([&minEdgeWeight = minEdgeWeight](node, node, Weight ew) {
                minEdgeWeight = std::min(minEdgeWeight, ew);
            });
        }

        if (G.isDirected()) {
            computeReachableNodesDir();
            G.forNodes([&](const auto u) { margGain[u] = closenessUBDir(u); });
        } else {
            pred.resize(n);
            computeReachableNodesUndir();
            G.forNodes([&](const auto u) { margGain[u] = closenessUBUnd(u); });
        }
    }

    void run() override;

    const std::vector<node> &groupMaxCloseness() const {
        assureFinished();
        return group;
    }

    // For debugging
    void checkDistFromGroup() const;

    std::vector<node> runNaive() {
        std::vector<node> result;
        std::fill(distFromGroup.begin(), distFromGroup.end(), infDist);

        while (result.size() < k) {
            node bestNode = none;
            double bestScore = -std::numeric_limits<double>::std::max();
            G.parallelForNodes([&](const auto u) {
                double curScore = 0;
                if (distFromGroup[u] != infDist)
                    curScore = -1. / (double)distFromGroup[u];

                const auto exploreNode = [&](node v, Weight d) -> void {
                    if (v == u)
                        return;
                    if (d < distFromGroup[v]) {
                        curScore += 1. / (double)d;
                        if (distFromGroup[v] != infDist)
                            curScore -= 1. / (double)distFromGroup[v];
                    }
#pragma omp critical
                    {
                        if (curScore > bestScore) {
                            bestScore = curScore;
                            bestNode = u;
                        }
                    }
                };

                if (G.isWeighted())
                    Traversal::DijkstraFrom(G, u, exploreNode);
                else
                    Traversal::BFSfrom(G, u, exploreNode);
            });

            assert(bestNode != none);
            result.push_back(bestNode);

            const auto updateDistances = [&distFromGroup = distFromGroup](
                                             node v, Weight d) -> void { distFromGroup[v] = d; };

            if (G.isWeighted())
                Traversal::DijkstraFrom(G, result.begin(), result.end(), updateDistances);
            else
                Traversal::BFSfrom(G, result.begin(), result.end(), updateDistances);
        }

        return result;
    }

    std::pair<bool, double> prunedDijkstraSwap(node source, double scoreDecr) {
        incTS();
        const uint8_t ts = tsGlobal[omp_get_thread_num()];
        auto &visited = visitedGlobal[omp_get_thread_num()];
        visited[source] = ts;
        auto &curDist = curDistGlobal[omp_get_thread_num()];
        std::fill(curDist.begin(), curDist.end(), infDist);
        curDist[source] = 0;
        auto &dijkstraHeap = dijkstraHeaps[omp_get_thread_num()];
        dijkstraHeap.clear();
        nearerNodes.clear();
        nearerNodes.push_back(source);

        const count reachableNodes = G.isDirected()
                                         ? reachableNodesUB[source]
                                         : reachableNodesInComponent[componentOfNode[source]];
        count visitedNodes = 1;

        const auto exploreNeighbors = [&](node u) -> void {
            for (const auto [w, weight] : G.weightNeighborRange(u)) {
                const edgeweight newDist = curDist[u] + weight;
                if (visited[w] != ts) {
                    visited[w] = ts;
                    if (newDist < nextDistFromGroup[w]) {
                        curDist[w] = newDist;
                        dijkstraHeap.push(w);
                    }
                } else if (newDist < std::min(nextDistFromGroup[w], curDist[w])) {
                    curDist[w] = newDist;
                    dijkstraHeap.update(w);
                }
            }
        };

        // Visit source now to avoid 'if (u != source)' in main loop
        exploreNeighbors(source);

        assert(nextDistFromGroup[source] > 0);
        double scoreIncr = -1. / nextDistFromGroup[source];

        group.push_back(source);
        assert(group.size() == k);
        std::vector<double> checkDist(G.upperNodeIdBound());
        for (node u : group)
            checkDist[u] = 0;
        Traversal::DijkstraFrom(G, group.begin(), group.end(),
                                [&](node u, edgeweight ew) { checkDist[u] = ew; });
        group.pop_back();

        while (!dijkstraHeap.empty()) {
            const node top = dijkstraHeap.extract_top();

            assert(curDist[top] == checkDist[top]);
            assert(nextDistFromGroup[top] > curDist[top]);
            assert(curDist[top] > 0);

            scoreIncr += 1. / curDist[top];
            if (nextDistFromGroup[top] != infDist)
                scoreIncr -= 1. / nextDistFromGroup[top];

            ++visitedNodes;
            const double scoreIncrUB =
                scoreIncr + static_cast<double>(reachableNodes - visitedNodes) / curDist[top];
            if (scoreIncrUB <= scoreDecr)
                return {false, 0};
            nearerNodes.push_back(top);

            exploreNeighbors(top);
        }

        return {true, scoreIncr};
    }

    std::pair<bool, double> prunedBFSSwap(node source, double scoreDecr) {
        incTS();
        auto &visited = visitedGlobal[omp_get_thread_num()];
        const uint8_t ts = tsGlobal[omp_get_thread_num()];
        visited[source] = ts;

        std::queue<node> q1, q2;
        q1.push(source);
        nearerNodes.clear();

        double scoreIncr = 0;
        if (nextDistFromGroup[source] != infDist)
            scoreIncr = -1. / static_cast<double>(nextDistFromGroup[source]);

        count level = 1;

        const count reachableNodes = G.isDirected()
                                         ? reachableNodesUB[source]
                                         : reachableNodesInComponent[componentOfNode[source]];
        count visitedNodes = 1;

        auto &curDist = curDistGlobal[omp_get_thread_num()];
        std::fill(curDist.begin(), curDist.end(), infDist);
        curDist[source] = 0;
        do {
            do {
                const node u = q1.front();
                q1.pop();
                nearerNodes.push_back(u);

                for (const node w : G.neighborRange(u)) {
                    if (visited[w] != ts) {
                        visited[w] = ts;
                        if (nextDistFromGroup[w] > level) {
                            curDist[w] = level;
                            q2.push(w);
                            scoreIncr += 1. / static_cast<double>(level);
                            if (nextDistFromGroup[w] != infDist)
                                scoreIncr -= 1. / static_cast<double>(nextDistFromGroup[w]);

                            ++visitedNodes;
                        }
                    }
                }
            } while (!q1.empty());

            ++level;
            const double scoreIncrUB = scoreIncr
                                       + static_cast<double>(reachableNodes - visitedNodes)
                                             / static_cast<double>(level + 1);
            if (scoreIncrUB <= scoreDecr)
                return {false, 0};
            std::swap(q1, q2);
        } while (!q1.empty());

        return {true, scoreIncr};
    }

    double computeScore() const {
        double score = 0;
        auto updateScore = [&](node, Weight w) -> void {
            if (w > 0)
                score += 1. / static_cast<double>(w);
        };
        if (G.isWeighted())
            Traversal::DijkstraFrom(G, group.begin(), group.end(), updateScore);
        else
            Traversal::BFSfrom(G, group.begin(), group.end(), updateScore);
        return score;
    }

    void checkNextDist() const {
        auto check = [&](node u, Weight ew) -> void {
            if (nextDistFromGroup[u] != ew)
                INFO(u, ": ", nextDistFromGroup[u], " -- ", ew);
            assert(nextDistFromGroup[u] == ew);
        };
        if (G.isWeighted())
            Traversal::DijkstraFrom(G, group.begin(), group.end(), check);
        else
            Traversal::BFSfrom(G, group.begin(), group.end(), check);
    }

    std::vector<node> runNaiveLs() {
        run();
        nextDistFromGroup.resize(G.upperNodeIdBound());
        std::vector<bool> inGroup(G.upperNodeIdBound());
        for (node u : group)
            inGroup[u] = true;

        double scoreDecr, scoreIncr;

        auto updateDistAfterRemoval = [&]() -> void {
            std::fill(nextDistFromGroup.begin(), nextDistFromGroup.end(), infDist);
            auto update = [&](node u, Weight w) -> void { nextDistFromGroup[u] = w; };
            assert(group.size() == k - 1);

            if (G.isWeighted())
                Traversal::DijkstraFrom(G, group.begin(), group.end(), update);
            else
                Traversal::BFSfrom(G, group.begin(), group.end(), update);

            scoreDecr = 0;
            G.parallelForNodes([&](node u) {
                if (distFromGroup[u] != nextDistFromGroup[u]) {
                    assert(nextDistFromGroup[u] != 0);
                    if (distFromGroup[u] > 0)
#pragma omp atomic update
                        scoreDecr += 1. / static_cast<double>(distFromGroup[u]);
                    if (nextDistFromGroup[u] != infDist)
#pragma omp atomic update
                        scoreDecr -= 1. / static_cast<double>(nextDistFromGroup[u]);
                }
            });
        };

        auto evalInsertion = [&](node v) -> bool {
            assert(!inGroup[v] && nextDistFromGroup[v] > 0);
            bool fullExploration;
            std::tie(fullExploration, scoreIncr) =
                G.isWeighted() ? prunedDijkstraSwap(v, scoreDecr) : prunedBFSSwap(v, scoreDecr);
#ifndef NDEBUG
            if (fullExploration) {
                double checkIncr = -computeScore();
                group.push_back(v);
                checkIncr += computeScore();
                group.pop_back();
                assert(std::abs(checkIncr - scoreIncr) <= 1e-6);
            }
#endif

            return (fullExploration && scoreIncr -scoreDecr > 1e-9);
        };

        bool anyProgress;
        do {
            anyProgress = false;
            node lastRemoved = none;

            assert(group.size() == k);
            for (size_t i = 0; i < group.size(); ++i) {
                lastRemoved = group[i];
                std::swap(group[i], group.back());
                group.pop_back();
                updateDistAfterRemoval();

#ifndef NDEBUG
                group.push_back(lastRemoved);
                double checkDecr = computeScore();
                group.pop_back();
                checkDecr -= computeScore();
                assert(std::abs(scoreDecr - checkDecr) <= 1e-6);
                checkNextDist();
#endif

                for (const node v : G.nodeRange()) {
                    if (inGroup[v])
                        continue;
                    if (evalInsertion(v)) {
                        anyProgress = true;
                        inGroup[lastRemoved] = false;
                        inGroup[v] = true;
                        group.push_back(v);
                        const auto &curDist = curDistGlobal[omp_get_thread_num()];
                        for (node x : nearerNodes) {
                            assert(curDist[x] < nextDistFromGroup[x]);
                            nextDistFromGroup[x] = curDist[x];
                        }
#ifndef NDEBUG
                        assert(group.size() == k);
                        assert(std::unordered_set<node>(group.begin(), group.end()).size() == k);
                        for (node u : group)
                            assert(inGroup[u]);
                        assert(std::count(inGroup.begin(), inGroup.end(), true) == k);
                        checkNextDist();
#endif
                        break;
                    }
                }

                if (!anyProgress) {
                    group.push_back(lastRemoved);
                    std::swap(group[i], group.back());
                } else {
                    std::swap(nextDistFromGroup, distFromGroup);
                    break;
                }
            }
        } while (anyProgress);

        assert(group.size() == k);
        return group;
    }
};

template <>
double GroupHarmonicClosenessGeneral<count>::closenessUBUnd(node u) const noexcept {
    const count reachableFromU = reachableNodesInComponent[componentOfNode[u]];
    const count degU = G.degree(u);
    double result = std::min(degU, reachableFromU);
    if (reachableFromU > degU + 1)
        result += static_cast<double>(reachableFromU - degU - 1) / 2.;
    return result;
}

template <>
double GroupHarmonicClosenessGeneral<edgeweight>::closenessUBUnd(node u) const noexcept {
    const count reachableFromU = reachableNodesInComponent[componentOfNode[u]];
    if (reachableFromU <= 1)
        return 0;
    edgeweight smallestWeight = infDist;
    G.forNeighborsOf(u, [&](node v, edgeweight ew) {
        if (ew < std::min(distFromGroup[v], smallestWeight))
            smallestWeight = ew;
    });
    return 1. / smallestWeight
           + static_cast<double>(reachableFromU - 2) / (smallestWeight + minEdgeWeight);
}

template <>
double GroupHarmonicClosenessGeneral<count>::closenessUBDir(node u) const noexcept {
    return G.degree(u) + static_cast<double>(reachableNodesUB[u] - G.degree(u) - 1) / 2.;
}

template <>
double GroupHarmonicClosenessGeneral<edgeweight>::closenessUBDir(node u) const noexcept {
    const count reachableFromU = reachableNodesUB[u];
    if (reachableFromU <= 1)
        return 0;
    const edgeweight smallestWeight =
        (*std::min_element(
             G.weightNeighborRange(u).begin(), G.weightNeighborRange(u).end(),
             [](const auto &e1, const auto &e2) -> bool { return e1.second < e2.second; }))
            .second;
    return static_cast<double>(G.degree(u)) / smallestWeight
           + static_cast<double>(reachableFromU - G.degree(u) - 1)
                 / (smallestWeight + minEdgeWeight);
}

template <class Weight>
void GroupHarmonicClosenessGeneral<Weight>::computeReachableNodesUndir() {
    ConnectedComponents cc(G);
    cc.run();

    reachableNodesInComponent.reserve(cc.numberOfComponents());
    for (const auto [cmpIndex, cmpSize] : cc.getComponentSizes())
        reachableNodesInComponent.push_back(cmpSize);
    G.forNodes([&](const auto u) { componentOfNode[u] = cc.componentOfNode(u); });
}

template <class Weight>
void GroupHarmonicClosenessGeneral<Weight>::computeReachableNodesDir() {
    reachableNodesUB.resize(G.upperNodeIdBound());
    StronglyConnectedComponents sccs(G);
    sccs.run();

    const count N = sccs.numberOfComponents();

    std::vector<count> reachL_scc(N, 0);
    std::vector<count> reachU_scc(N, 0);
    std::vector<count> reachU_without_max_scc(N, 0);
    std::vector<bool> reach_from_max_scc(N, false);
    std::vector<bool> reaches_max_scc(N, false);
    std::vector<std::vector<count>> sccs_vec(N, std::vector<count>());
    Graph sccGraph(N, false, true);
    std::vector<bool> found(N, false);
    count maxSizeCC = 0;

    // We compute the vector sccs_vec, where each component contains the list of
    // its nodes
    G.forNodes([&](const auto v) { sccs_vec[sccs.componentOfNode(v)].push_back(v); });

    // We compute the SCC graph and store it in sccGraph
    for (count V = 0; V < N; V++) {
        for (count v : sccs_vec[V]) {
            G.forNeighborsOf(v, [&](node w) {
                count W = sccs.componentOfNode(w);

                if (W != V && !found[W]) {
                    found[W] = true;
                    sccGraph.addEdge(V, W);
                }
            });
        }
        sccGraph.forNeighborsOf(V, [&](node W) { found[W] = false; });
        if (sccGraph.degreeOut(V) > sccGraph.degreeOut(maxSizeCC))
            maxSizeCC = V;
    }

    // BFS from the biggest SCC.
    std::queue<count> Q;
    Q.push(maxSizeCC);
    reach_from_max_scc[maxSizeCC] = true;
    do {
        count V = Q.front();
        Q.pop();
        reachL_scc[maxSizeCC] += sccs_vec[V].size();
        sccGraph.forNeighborsOf(V, [&](node W) {
            if (!reach_from_max_scc[W]) {
                reach_from_max_scc[W] = true;
                Q.push(W);
            }
        });
    } while (!Q.empty());

    reachU_scc[maxSizeCC] = reachL_scc[maxSizeCC];
    reaches_max_scc[maxSizeCC] = true;

    // so far only the largest SCC has reach_U and reach_L > 0

    // Dynamic programming to compute number of reachable vertices
    for (count V = 0; V < N; V++) {
        if (V == maxSizeCC)
            continue;

        sccGraph.forNeighborsOf(V, [&](node W) {
            reachL_scc[V] = std::max(reachL_scc[V], reachL_scc[W]);
            if (!reach_from_max_scc[W])
                reachU_without_max_scc[V] += reachU_without_max_scc[W];

            reachU_scc[V] += reachU_scc[W];
            reachU_scc[V] = std::min(reachU_scc[V], G.numberOfNodes());
            reaches_max_scc[V] = reaches_max_scc[V] || reaches_max_scc[W];
        });

        if (reaches_max_scc[V])
            reachU_scc[V] = reachU_without_max_scc[V] + reachU_scc[V];

        reachL_scc[V] += sccs_vec[V].size();
        reachU_scc[V] += sccs_vec[V].size();
        reachU_scc[V] = std::min(reachU_scc[V], G.numberOfNodes());
    }

    G.forNodes([&](const auto v) { reachableNodesUB[v] = reachU_scc[sccs.componentOfNode(v)]; });
}

template <>
std::pair<bool, double> GroupHarmonicClosenessGeneral<count>::bfscut(node source,
                                                                     double bestScore) {
    incTS();
    const uint8_t ts = tsGlobal[omp_get_thread_num()];
    auto &visited = visitedGlobal[omp_get_thread_num()];
    visited[source] = ts;
    std::queue<node> q1, q2;
    q1.push(source);
    const count reachableFromSource = G.isDirected()
                                          ? reachableNodesUB[source]
                                          : reachableNodesInComponent[componentOfNode[source]];
    double curScore = 0, scoreUB = margGain[source];
    count level = 1, visitedNodes = 1;
    const count undirected = !G.isDirected();

    auto &curDist = curDistGlobal[omp_get_thread_num()];
    std::fill(curDist.begin(), curDist.end(), infDist);
    curDist[source] = 0;
    do {
        count nodesAtNextLevelUB = 0, prevAtNextLevelUB = 0;
        do {
            const node u = q1.front();
            q1.pop();

            for (const node w : G.neighborRange(u)) {
                if (visited[w] != ts) {
                    visited[w] = ts;
                    curDist[w] = curDist[u] + 1;
                    q2.push(w);
                    ++visitedNodes;
                    curScore += 1. / static_cast<double>(level);
                    nodesAtNextLevelUB += G.degree(w) - undirected;
                    if (undirected)
                        pred[w] = u;
                } else if (undirected && prevAtNextLevelUB && u != pred[w]) {
                    --prevAtNextLevelUB;
                    assert(u != source);
                    scoreUB -=
                        1. / static_cast<double>(level) - 1. / static_cast<double>(level + 1);
                    assert(scoreUB > 0);
                    if (scoreUB <= bestScore)
                        return {false, scoreUB};
                }
            }
        } while (!q1.empty());

        assert(visitedNodes <= reachableFromSource);
        nodesAtNextLevelUB = std::min(nodesAtNextLevelUB, reachableFromSource - visitedNodes);
        prevAtNextLevelUB = nodesAtNextLevelUB;
        scoreUB = curScore
                  + static_cast<double>(nodesAtNextLevelUB) / static_cast<double>(level + 1)
                  + static_cast<double>(reachableFromSource - visitedNodes - nodesAtNextLevelUB)
                        / static_cast<double>(level + 2);

        if (scoreUB <= bestScore)
            return {false, scoreUB};

        ++level;
        std::swap(q2, q1);
    } while (!q1.empty());

    return {true, curScore};
}

template <>
std::pair<bool, double> GroupHarmonicClosenessGeneral<edgeweight>::bfscut(node source,
                                                                          double bestScore) {
    incTS();
    auto &visited = visitedGlobal[omp_get_thread_num()];
    const uint8_t ts = tsGlobal[omp_get_thread_num()];
    visited[source] = ts;

    const count reachableFromSource = G.isDirected()
                                          ? reachableNodesUB[source]
                                          : reachableNodesInComponent[componentOfNode[source]];
    double curScore = 0, scoreUB = 0;
    count visitedNodes = 1;

    auto &dijkstraHeap = dijkstraHeaps[omp_get_thread_num()];
    auto &curDist = curDistGlobal[omp_get_thread_num()];

    const auto exploreNeighbors = [&](node u) -> void {
        for (const auto [w, weight] : G.weightNeighborRange(u)) {
            const edgeweight newDist = curDist[u] + weight;
            if (visited[w] != ts) {
                visited[w] = ts;
                curDist[w] = newDist;
                dijkstraHeap.push(w);
            } else if (newDist < curDist[w]) {
                curDist[w] = newDist;
                dijkstraHeap.update(w);
            }
        }
    };

    std::fill(curDist.begin(), curDist.end(), infDist);
    curDist[source] = 0;
    dijkstraHeap.clear();
    // Explore source now to avoid "if (u != source) curScore += 1. / curDist[u];" in main
    // loop
    exploreNeighbors(source);
    do {
        const node u = dijkstraHeap.extract_top();
        curScore += 1. / curDist[u];
        ++visitedNodes;
        scoreUB = curScore
                  + static_cast<double>(reachableFromSource - visitedNodes)
                        / (curDist[u] + minEdgeWeight);
        if (scoreUB <= bestScore)
            return {false, scoreUB};
        exploreNeighbors(u);
    } while (!dijkstraHeap.empty());

    return {true, curScore};
}

template <>
std::pair<bool, double> GroupHarmonicClosenessGeneral<count>::prunedBFS(node source,
                                                                        double bestScore) {
    incTS();
    auto &visited = visitedGlobal[omp_get_thread_num()];
    const uint8_t ts = tsGlobal[omp_get_thread_num()];
    visited[source] = ts;

    std::queue<node> q1, q2;
    q1.push(source);
    auto &curNearerNodes = curNearerNodesGlobal[omp_get_thread_num()];
    curNearerNodes.clear();

    const count reachableFromSource = G.isDirected()
                                          ? reachableNodesUB[source]
                                          : reachableNodesInComponent[componentOfNode[source]];
    double curScore = 0, scoreUB = margGain[source];
    if (distFromGroup[source] != infDist) {
        curScore = -1. / static_cast<double>(distFromGroup[source]);
        scoreUB -= 1. / static_cast<double>(distFromGroup[source]);
    }

    count level = 1, visitedNodes = (distFromGroup[source] != 1);
    const count undirected = !G.isDirected();

    auto &curDist = curDistGlobal[omp_get_thread_num()];
    std::fill(curDist.begin(), curDist.end(), infDist);
    curDist[source] = 0;
    do {
        count nodesAtNextLevelUB = 0, prevAtNextLevelUB = 0;
        do {
            const node u = q1.front();
            q1.pop();
            curNearerNodes.push_back(u);

            for (const node w : G.neighborRange(u)) {
                if (visited[w] != ts) {
                    visited[w] = ts;
                    ++visitedNodes;
                    if (undirected)
                        pred[w] = u;
                    if (distFromGroup[w] > level) {
                        curDist[w] = level;
                        q2.push(w);
                        nodesAtNextLevelUB += G.degree(w) - undirected;
                        curScore += 1. / static_cast<double>(level);
                        if (distFromGroup[w] != infDist)
                            curScore -= 1. / static_cast<double>(distFromGroup[w]);
                    }
                } else if (undirected && prevAtNextLevelUB && u != pred[w]) {
                    --prevAtNextLevelUB;
                    scoreUB -=
                        1. / static_cast<double>(level) - 1. / static_cast<double>(level + 1);
                    assert(scoreUB > 0);
                    if (scoreUB <= bestScore)
                        return {false, scoreUB};
                }
            }
        } while (!q1.empty());

        scoreUB =
            curScore + static_cast<double>(nodesAtNextLevelUB) / static_cast<double>(level + 1);
        if (reachableFromSource > visitedNodes + nodesAtNextLevelUB)
            scoreUB += static_cast<double>(reachableFromSource - visitedNodes - nodesAtNextLevelUB)
                       / static_cast<double>(level + 2);

        if (scoreUB <= bestScore)
            return {false, scoreUB};

        ++level;
        std::swap(q1, q2);
    } while (!q1.empty());

    return {true, curScore};
}

template <>
std::pair<bool, double> GroupHarmonicClosenessGeneral<edgeweight>::prunedBFS(node source,
                                                                             double bestScore) {
    incTS();
    auto &visited = visitedGlobal[omp_get_thread_num()];
    const uint8_t ts = tsGlobal[omp_get_thread_num()];
    auto &dijkstraHeap = dijkstraHeaps[omp_get_thread_num()];
    dijkstraHeap.clear();
    visited[source] = ts;
    auto &curNearerNodes = curNearerNodesGlobal[omp_get_thread_num()];
    curNearerNodes.clear();
    curNearerNodes.push_back(source);

    const count reachableFromSource = G.isDirected()
                                          ? reachableNodesUB[source]
                                          : reachableNodesInComponent[componentOfNode[source]];
    double curScore = 0, scoreUB = margGain[source];
    assert(distFromGroup[source] > 0);
    if (distFromGroup[source] != infDist) {
        curScore = -1. / distFromGroup[source];
        scoreUB -= 1. / distFromGroup[source];
    }

    count visitedNodes = 1;
    auto &curDist = curDistGlobal[omp_get_thread_num()];
    std::fill(curDist.begin(), curDist.end(), infDist);
    curDist[source] = 0.;

    const auto exploreNeighbors = [&](node u) -> void {
        for (const auto [w, weight] : G.weightNeighborRange(u)) {
            const edgeweight newDist = curDist[u] + weight;
            if (visited[w] != ts) {
                visited[w] = ts;
                if (newDist < distFromGroup[w]) {
                    curDist[w] = newDist;
                    dijkstraHeap.push(w);
                }
            } else if (newDist < std::min(distFromGroup[w], curDist[w])) {
                curDist[w] = newDist;
                dijkstraHeap.update(w);
            }
        }
    };

    // Visit source now to avoid 'if (u != source)' in main loop
    exploreNeighbors(source);

    if (dijkstraHeap.empty())
        return {true, -1. / distFromGroup[source]};

    do {
        const node u = dijkstraHeap.extract_top();
        ++visitedNodes;
        curNearerNodes.push_back(u);
        assert(curDist[u] < distFromGroup[u]);
        curScore += 1. / curDist[u];
        if (distFromGroup[u] != infDist)
            curScore -= 1. / distFromGroup[u];

        scoreUB = curScore + static_cast<double>(reachableFromSource - visitedNodes) / curDist[u];

        if (scoreUB <= bestScore)
            return {false, scoreUB};
        exploreNeighbors(u);
    } while (!dijkstraHeap.empty());

    return {true, curScore};
}

template <class Weight>
node GroupHarmonicClosenessGeneral<Weight>::findTopHCloseness() {
    node bestNode = none;
    double bestScore = 0;
    candidateNodesPQ.build_heap(G.nodeRange().begin(), G.nodeRange().end());

    std::atomic<bool> stop{false};

#pragma omp parallel
    {
        while (!stop.load(std::memory_order_relaxed)) {
            node u = none;

#pragma omp critical
            {
                if (candidateNodesPQ.empty()) {
                    stop.store(true, std::memory_order_relaxed);
                } else {
                    u = candidateNodesPQ.extract_top();
                    if (margGain[u] <= bestScore) {
                        stop.store(true, std::memory_order_relaxed);
                        u = none;
                    }
                }
            }

            if (u == none)
                break;

            bool exact;
            double score;
            std::tie(exact, score) = bfscut(u, bestScore);
            margGain[u] = score;
#pragma omp critical
            {
                if (exact && score > bestScore) {
                    bestNode = u;
                    bestScore = score;
                    std::swap(curDistGlobal[omp_get_thread_num()], distFromGroup);
#ifndef NDEBUG
                    group.push_back(bestNode);
                    checkDistFromGroup();
                    group.clear();
#endif
                }
            }
        }
    }

    if (!G.isDirected() && !G.isDirected())
        reachableNodesInComponent[componentOfNode[bestNode]] -= G.degree(bestNode) + 1;
    return bestNode;
}

template <class Weight>
node GroupHarmonicClosenessGeneral<Weight>::findNodeWithHighestMargGain() {
    node bestNode = none;
    double bestScore = -std::numeric_limits<double>::std::max();
    if (!G.isDirected())
        G.forNodes([&](const auto u) {
            if (distFromGroup[u] > 0)
                margGain[u] = std::min(margGain[u], closenessUBUnd(u));
        });

    candidateNodesPQ.build_heap(G.nodeRange().begin(), G.nodeRange().end());
    for (const node u : group)
        candidateNodesPQ.remove(u);

    std::atomic<bool> stop{false};
#pragma omp parallel
    {
        while (!stop.load(std::memory_order_relaxed)) {
            node u = none;
#pragma omp critical
            {
                if (candidateNodesPQ.empty()) {
                    stop.store(true, std::memory_order_relaxed);
                } else {
                    u = candidateNodesPQ.extract_top();
                    if (margGain[u] <= bestScore) {
                        stop.store(true, std::memory_order_relaxed);
                        u = none;
                    }
                }
            }

            if (u == none)
                break;

            bool exact;
            double score;
            std::tie(exact, score) = prunedBFS(u, bestScore);
            margGain[u] = score;

#pragma omp critical
            {
                if (exact && score > bestScore) {
                    bestNode = u;
                    bestScore = score;
                    std::swap(curDistGlobal[omp_get_thread_num()], curBestDist);
                    std::swap(curNearerNodesGlobal[omp_get_thread_num()], nearerNodes);
                }
            }
        }
    }

    assert(bestNode != none);
    for (const node u : nearerNodes) {
        assert(distFromGroup[u] > curBestDist[u]);
        distFromGroup[u] = curBestDist[u];
        if (!G.isWeighted() && !G.isDirected() && distFromGroup[u] > 1 && curBestDist[u] <= 1) {
            assert(reachableNodesInComponent[componentOfNode[u]]);
            --reachableNodesInComponent[componentOfNode[u]];
        }
    }

    return bestNode;
}

template <class Weight>
void GroupHarmonicClosenessGeneral<Weight>::run() {
    group.push_back(findTopHCloseness());
#ifndef NDEBUG
    checkDistFromGroup();
#endif

    while (group.size() < k) {
        group.push_back(findNodeWithHighestMargGain());
#ifndef NDEBUG
        checkDistFromGroup();
#endif
        margGain[group.back()] = 0;
    }

    hasRun = true;
}

template <>
void GroupHarmonicClosenessGeneral<count>::checkDistFromGroup() const {
    Traversal::BFSfrom(G, group.begin(), group.end(), [&](node x, count d) {
        if (distFromGroup[x] != d)
            INFO("Error in node ", x, ": ", distFromGroup[x], " -- ", d);
        assert(distFromGroup[x] == d);
    });
}

template <>
void GroupHarmonicClosenessGeneral<edgeweight>::checkDistFromGroup() const {
    Traversal::DijkstraFrom(G, group.begin(), group.end(), [&](node x, edgeweight d) {
        if (distFromGroup[x] != d)
            INFO("Error in node ", x, ": ", distFromGroup[x], " -- ", d);
        assert(distFromGroup[x] == d);
    });
}

using GroupHarmonicCloseness = GroupHarmonicClosenessGeneral<count>;
using GroupHarmonicClosenessWeighted = GroupHarmonicClosenessGeneral<edgeweight>;

} // namespace NetworKit
#endif

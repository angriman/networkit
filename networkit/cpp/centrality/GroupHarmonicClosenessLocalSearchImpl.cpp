/*
 *  GroupHarmonicClosenessLocalSearchImpl.cpp
 *
 *  Created on: 20.12.2020
 *     Authors: Eugenio Angriman <angrimae@hu-berlin.de>
 *              Alexander van der Grinten <avdgrinten@hu-berlin.de>
 */

// networkit-format

#include <cassert>
#include <omp.h>
#include <stdexcept>

#include <networkit/graph/BFS.hpp>
#include <networkit/graph/Dijkstra.hpp>
#include <networkit/graph/GraphTools.hpp>

#include "GroupHarmonicClosenessLocalSearchImpl.hpp"

namespace NetworKit {
namespace GroupHarmonicClosenessLocalSearchDetails {

template <class Weight>
GroupHarmonicClosenessLocalSearchImpl<Weight>::GroupHarmonicClosenessLocalSearchImpl(const Graph &G)
    : G(&G) {

    const count n = G.upperNodeIdBound(), threads = omp_get_max_threads();

    inGroup.resize(n, false);
    visitedGlobal.resize(threads, std::vector<bool>(n, false));
    distFromGroup.resize(n);
    newDistFromGroup.resize(n);
    distGlobal.resize(threads, std::vector<Weight>(n, infDist));

    if (G.isWeighted()) {
        dijkstraHeaps.reserve(threads);
        for (int i = 0; i < omp_get_max_threads(); ++i)
            dijkstraHeaps.emplace_back(Less{distGlobal[i]});
    }
}

template <class Weight>
GroupHarmonicClosenessLocalSearchImpl<Weight>::GroupHarmonicClosenessLocalSearchImpl(const Graph &G,
                                                                                     count k)
    : GroupHarmonicClosenessLocalSearchImpl(G), k(k) {

    if (k == 0 || k >= G.numberOfNodes())
        throw std::runtime_error("Error, k must be in [1, n - 1].");

    group = generateRandomGroup(k);
}

template <class Weight>
template <class InputIt>
GroupHarmonicClosenessLocalSearchImpl<Weight>::GroupHarmonicClosenessLocalSearchImpl(const Graph &G,
                                                                                     InputIt first,
                                                                                     InputIt last)
    : GroupHarmonicClosenessLocalSearchImpl(G), group(first, last), k(group.size()) {

    for (node u : group)
        inGroup[u] = true;
}

template <class Weight>
std::vector<node> GroupHarmonicClosenessLocalSearchImpl<Weight>::generateRandomGroup() {
    std::vector<node> randomGroup;
    randomGroup.reserve(k);

    do {
        const node randomNode = GraphTools::randomNode(*G);
        if (!inGroup[randomNode]) {
            randomGroup.push_back(randomNode);
            inGroup[randomNode] = true;
        }
    } while (randomGroup.size() < k);

    return randomGroup;
}

template <class Weight>
GroupHarmonicClosenessLocalSearchImpl<Weight>::~GroupHarmonicClosenessLocalSearchImpl() = default;

template <class Weight>
void GroupHarmonicClosenessLocalSearchImpl<Weight>::run(count maxSwaps) {
    INFO(maxSwaps);
    computeDistFromGroup();

    double scoreDecrement, scoreIncrement;

    const auto updateDistAfterRemoval = [&]() -> void {
        std::fill(newDistFromGroup.begin(), newDistFromGroup.end(), infDist);
        assert(group.size() == k - 1);

        computeDistFromGroup();

        scoreDecrement = 0;
        G.parallelForNodes([&](node u) {
            if (distFromGroup[u] != newDistFromGroup[u]) {
                assert(newDistFromGroup[u] != 0);
                if (distFromGroup[u] > 0)
#pragma omp atomic update
                    scoreDecrement += 1. / static_cast<double>(distFromGroup[u]);
                if (newDistFromGroup[u] != infDist)
#pragma omp atomic update
                    scoreDecrement -= 1. / static_cast<double>(newDistFromGroup[u]);
            }
        });
    };

        const auto evalInsertion = [&](node v) -> bool {
            assert(!inGroup[v] && newDistFromGroup[v] > 0);
            auto ssspResult = prunedSSSP(v, scoreDecr);
#ifndef NDEBUG
            if (!ssspResult.pruned) {
                double checkIncr = -computeScore();
                group.push_back(v);
                checkIncr += computeScore();
                group.pop_back();
                assert(std::abs(checkIncr - scoreIncr) <= 1e-6);
            }
#endif

            return fullExploration && scoreIncr - scoreDecr > 0;
        };
}

template <>
void GroupHarmonicClosenessLocalSearchImpl<count>::computeDistFromGroup() {
    Traversal::BFSfrom(
        *G, group.begin(), group.end(),
        [&distFromGroup = distFromGroup](node u, count dist) { distFromGroup[u] = dist; });
}

template <>
void GroupHarmonicClosenessLocalSearchImpl<edgeweight>::computeDistFromGroup() {
    Traversal::DijkstraFrom(
        *G, group.begin(), group.end(),
        [&distFromGroup = distFromGroup](node u, edgeweight dist) { distFromGroup[u] = dist; });
}

template class GroupHarmonicClosenessLocalSearchImpl<count>;
template class GroupHarmonicClosenessLocalSearchImpl<edgeweight>;

} // namespace GroupHarmonicClosenessLocalSearchDetails
} // namespace NetworKit

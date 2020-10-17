/*
 *  MultiSourceDijkstra.cpp
 *  Created on: 17.10.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#include <networkit/distance/MultiSourceDijkstra.hpp>

namespace NetworKit {

void MultiSourceDijkstra::run() {
    dist.resize(G->upperNodeIdBound());
    std::fill(dist.begin(), dist.end(), infdist);
    reachedNodes = 0;

    for (node source : sources)
        dist[source] = 0;

    heap.clear();
    heap.reserve(G->upperNodeIdBound());
    heap.build_heap(sources.begin(), sources.end());

    if (targets.empty())
        runWithoutTargets();
    else
        runWithTargets();

    hasRun = true;
}

void MultiSourceDijkstra::runWithoutTargets() {
    while (!heap.empty()) {
        const node u = heap.extract_top();
        ++reachedNodes;

        G->forNeighborsOf(u, [&](node v, edgeweight ew) {
            const edgeweight newDist = dist[u] + ew;
            if (newDist < dist[v]) {
                dist[v] = newDist;
                heap.update(v);
            }
        });
    }
}

void MultiSourceDijkstra::runWithTargets() {
    std::unordered_set<node> targetsToVisit;
    for (node target : targets)
        if (dist[target] == infdist)
            targetsToVisit.insert(target);

    while (!heap.empty()) {
        const node u = heap.extract_top();

        auto it = targetsToVisit.find(u);
        if (it != targetsToVisit.end()) {
            targetsToVisit.erase(it);
            if (targetsToVisit.empty())
                break;
        }

        G->forNeighborsOf(u, [&](node v, edgeweight ew) {
            const edgeweight newDist = dist[u] + ew;
            if (newDist < dist[v]) {
                dist[v] = newDist;
                heap.update(v);
                reachedNodes += static_cast<count>(dist[v] == infdist);
            }
        });
    }
}

} // namespace NetworKit

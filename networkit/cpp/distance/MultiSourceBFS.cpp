/*
 *  MultiSourceBFS.cpp
 *  Created on: 17.10.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#include <algorithm>
#include <queue>

#include <networkit/distance/MultiSourceBFS.hpp>

namespace NetworKit {

constexpr edgeweight MSSP::infdist;

void MultiSourceBFS::run() {
    dist.resize(G->upperNodeIdBound());
    std::fill(dist.begin(), dist.end(), infdist);
    reachedNodes = 0;

    std::queue<node> queue;
    for (node u : sources) {
        queue.push(u);
        dist[u] = 0;
    }

    auto runWithoutTargets = [&]() -> void {
        while (!queue.empty()) {
            const node u = queue.front();
            queue.pop();
            ++reachedNodes;

            G->forNeighborsOf(u, [&](node v) {
                if (dist[v] == infdist) {
                    queue.push(v);
                    dist[v] = dist[u] + 1;
                }
            });
        }
    };

    auto runWithTargets = [&]() -> void {
        std::unordered_set<node> targetsToVisit;
        for (node target : targets)
            if (dist[target] == infdist)
                targetsToVisit.insert(target);

        while (!queue.empty() && !targetsToVisit.empty()) {
            const node u = queue.front();
            queue.pop();

            for (node v : G->neighborRange(u)) {
                if (dist[v] == infdist) {
                    queue.push(v);
                    dist[v] = dist[u] + 1;
                    ++reachedNodes;

                    auto it = targetsToVisit.find(v);
                    if (it != targetsToVisit.end()) {
                        targetsToVisit.erase(it);
                        if (targetsToVisit.empty())
                            break;
                    }
                }
            }
        }
    };

    if (targets.empty())
        runWithoutTargets();
    else
        runWithTargets();

    hasRun = true;
}

} // namespace NetworKit

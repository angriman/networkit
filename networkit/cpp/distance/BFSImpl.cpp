#include "BFSDetails.hpp"

#include <queue>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace NetworKit::BFSDetails {

void BFSImpl::run() {
    const count z = G->upperNodeIdBound();
    reachedNodes = 1;
    sumDist = 0.;

    std::fill(distances.begin(), distances.end(), infDist);

    if (distances.size() < z)
        distances.resize(z, infDist);

    if (storePaths) {
        previous.clear();
        previous.resize(z);
        npaths.clear();
        npaths.resize(z, 0);
        npaths[source] = 1;
    }

    if (storeNodesSortedByDistance) {
        std::vector<node> empty;
        std::swap(nodesSortedByDistance, empty);
    }

    std::queue<node> q;
    q.push(source);
    distances[source] = 0.;

    const auto visitNeighbors = [&](node v) -> void {
        if (distances[v] == infDist) {
            q.push(v);
            distances[v] = distances[u] + 1.;
            sumDist += distances[v];
            ++reachedNodes;
            if (storePaths) {
                previous[v] = {u};
                npaths[v] = npaths[u];
            }
        } else if (storePaths && (distances[v] == distances[u] + 1.)) {
            // additional predecessor
            previous[v].push_back(u);
            // all the shortest paths to u are also  shortest paths to v now
            npaths[v] += npaths[u];
        }
    };

    while (!q.empty()) {
        node u = q.front();
        q.pop();

        if (storeNodesSortedByDistance)
            nodesSortedByDistance.push_back(u);

        if (target == u)
            break;

        // insert untouched neighbors into queue
        if (Reverse)
            G->forInNeighborsOf(u, visitNeighbors);
        else
            G->forNeighborsOf(u, visitNeighbors);
    }

    hasRun = true;
}

} // namespace NetworKit::BFSDetails

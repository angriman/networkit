/*
 *  MSSP.hpp
 *  Created on: 17.10.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/distance/MSSP.hpp>

namespace NetworKit {
const std::vector<edgeweight> &MSSP::getDistances() const {
    assureFinished();
    return dist;
}

edgeweight MSSP::distance(node target) const {
    assureFinished();
    return dist[target];
}

std::vector<node> MSSP::getNodesSortedByDistance() const {
    assureFinished();
    std::vector<node> nodesSortedByDistance(G->nodeRange().begin(), G->nodeRange().end());
    Aux::Parallel::sort(nodesSortedByDistance.begin(), nodesSortedByDistance.end(),
                        [&](node u, node v) { return dist[u] < dist[v]; });
    return nodesSortedByDistance;
}

std::vector<node> MSSP::getTargetNodesSortedByDistance() const {
    assureFinished();
    std::vector<node> sortedTargetNodes(targets);
    Aux::Parallel::sort(sortedTargetNodes.begin(), sortedTargetNodes.end(),
                        [&](node u, node v) { return dist[u] < dist[v]; });
    return sortedTargetNodes;
}

std::unordered_map<node, edgeweight> MSSP::getDistancesToTargets() const {
    assureFinished();
    std::unordered_map<node, edgeweight> distToTargets;
    for (node target : targets)
        distToTargets[target] = dist[target];

    return distToTargets;
}

count MSSP::getReachableNodes() const {
    assureFinished();
    return reachedNodes;
}

} // namespace NetworKit

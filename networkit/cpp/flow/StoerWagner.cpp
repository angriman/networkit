// networkit-format
#include <algorithm>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/flow/StoerWagner.hpp>

namespace NetworKit {

StoerWagner::StoerWagner(const Graph &G) : G(&G) {
    if (G.isDirected()) {
        throw std::runtime_error("StoerWagner: the graph must be undirected.");
    }
    if (G.numberOfNodes() < 2) {
        throw std::runtime_error("StoerWagner: the graph must have at least two nodes.");
    }
    G.forEdges([&](node, node, edgeweight weight) -> void {
        if (weight < 0) {
            throw std::runtime_error("StoerWagner: edge weights must not be negative.");
        }
    });
}

void StoerWagner::run() {
    Aux::SignalHandler handler;

    // Node ids may have gaps, super-nodes are consecutive.
    std::vector<index> superNodeOf(G->upperNodeIdBound(), none);
    std::vector<std::vector<node>> members;
    members.reserve(G->numberOfNodes());
    G->forNodes([&](node u) -> void {
        superNodeOf[u] = members.size();
        members.push_back({u});
    });

    std::vector<std::unordered_map<index, edgeweight>> neighbors(members.size());
    G->forEdges([&](node u, node v, edgeweight weight) -> void {
        if (u != v) {
            neighbors[superNodeOf[u]][superNodeOf[v]] += weight;
            neighbors[superNodeOf[v]][superNodeOf[u]] += weight;
        }
    });

    cutValue = std::numeric_limits<edgeweight>::infinity();
    cutSide.clear();

    std::vector<index> active(members.size());
    for (index i = 0; i < active.size(); ++i) {
        active[i] = i;
    }
    std::vector<edgeweight> attachment(members.size());
    std::vector<unsigned char> added(members.size());

    while (active.size() > 1) {
        handler.assureRunning();

        // Add the super-nodes one by one, always picking the most tightly connected one. Outdated
        // queue entries are skipped when popped.
        std::priority_queue<std::pair<edgeweight, index>> queue;
        for (const index superNode : active) {
            attachment[superNode] = 0;
            added[superNode] = 0;
            queue.emplace(0, superNode);
        }

        index previous = none, last = none;
        for (count numberAdded = 0; numberAdded < active.size();) {
            const auto [priority, current] = queue.top();
            queue.pop();
            if (added[current] || priority != attachment[current]) {
                continue;
            }
            added[current] = 1;
            ++numberAdded;
            previous = last;
            last = current;

            for (const auto &[neighbor, weight] : neighbors[current]) {
                if (!added[neighbor]) {
                    attachment[neighbor] += weight;
                    queue.emplace(attachment[neighbor], neighbor);
                }
            }
        }

        // The cut between the last super-node and everything else is the minimum cut between the
        // last two super-nodes, so they can be merged unless that is the overall minimum.
        if (attachment[last] < cutValue) {
            cutValue = attachment[last];
            cutSide = members[last];
        }

        if (neighbors[previous].size() < neighbors[last].size()) {
            std::swap(previous, last);
        }
        neighbors[previous].erase(last);
        for (const auto &[neighbor, weight] : neighbors[last]) {
            if (neighbor != previous) {
                neighbors[previous][neighbor] += weight;
                neighbors[neighbor][previous] += weight;
                neighbors[neighbor].erase(last);
            }
        }
        neighbors[last].clear();

        if (members[previous].size() < members[last].size()) {
            std::swap(members[previous], members[last]);
        }
        members[previous].insert(members[previous].end(), members[last].begin(),
                                 members[last].end());
        members[last].clear();

        active.erase(std::find(active.begin(), active.end(), last));
    }

    std::sort(cutSide.begin(), cutSide.end());
    hasRun = true;
}

} // namespace NetworKit

// networkit-format
#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/SteinerTree.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace NetworKit {

namespace {

// An edge between the regions of two terminals, and the length of the path through it.
struct Bridge {
    edgeweight length;
    node u;
    node v;
    edgeweight weight;
};

} // namespace

SteinerTree::SteinerTree(const Graph &G, std::vector<node> terminals)
    : G(&G), terminals(std::move(terminals)) {
    if (G.isDirected()) {
        throw std::runtime_error("SteinerTree: the graph must be undirected.");
    }
    G.forEdges([](node, node, edgeweight weight) -> void {
        if (weight < 0) {
            throw std::runtime_error("SteinerTree: edge weights must not be negative.");
        }
    });
    for (const node terminal : this->terminals) {
        if (!G.hasNode(terminal)) {
            throw std::runtime_error("SteinerTree: the terminals must be in the graph.");
        }
    }
    std::sort(this->terminals.begin(), this->terminals.end());
    this->terminals.erase(std::unique(this->terminals.begin(), this->terminals.end()),
                          this->terminals.end());
}

void SteinerTree::run() {
    Aux::SignalHandler handler;
    const count z = G->upperNodeIdBound();

    connected = false;
    totalWeight = 0;
    tree = Graph(z, true);
    for (node u = 0; u < z; ++u) {
        if (!G->hasNode(u)) {
            tree.removeNode(u);
        }
    }

    // Assign each node to its closest terminal. This is a Dijkstra search from all of them.
    std::vector<edgeweight> distance(z, std::numeric_limits<edgeweight>::infinity());
    std::vector<node> closestTerminal(z, none);
    std::vector<node> parent(z, none);
    std::vector<edgeweight> parentWeight(z, 0);

    using QueueEntry = std::pair<edgeweight, node>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<>> queue;
    for (const node terminal : terminals) {
        distance[terminal] = 0;
        closestTerminal[terminal] = terminal;
        queue.emplace(0., terminal);
    }

    while (!queue.empty()) {
        const edgeweight reached = queue.top().first;
        const node u = queue.top().second;
        queue.pop();
        if (reached > distance[u]) {
            continue;
        }
        G->forNeighborsOf(u, [&](node v, edgeweight weight) -> void {
            if (reached + weight < distance[v]) {
                distance[v] = reached + weight;
                closestTerminal[v] = closestTerminal[u];
                parent[v] = u;
                parentWeight[v] = weight;
                queue.emplace(distance[v], v);
            }
        });
    }
    handler.assureRunning();

    std::vector<Bridge> bridges;
    G->forEdges([&](node u, node v, edgeweight weight) -> void {
        if (closestTerminal[u] != none && closestTerminal[v] != none
            && closestTerminal[u] != closestTerminal[v]) {
            bridges.push_back({distance[u] + weight + distance[v], u, v, weight});
        }
    });
    std::sort(bridges.begin(), bridges.end(), [](const Bridge &a, const Bridge &b) -> bool {
        return std::tie(a.length, a.u, a.v) < std::tie(b.length, b.u, b.v);
    });

    // The bridges and the paths to the terminals cannot form cycles: paths stay within a region,
    // where they are part of a tree, and the selected bridges connect regions without cycles.
    std::vector<bool> connectedToTerminal(z, false);
    const auto addEdge = [&](node u, node v, edgeweight weight) -> void {
        tree.addEdge(u, v, weight);
        totalWeight += weight;
    };
    const auto addPathToTerminal = [&](node u) -> void {
        while (parent[u] != none && !connectedToTerminal[u]) {
            connectedToTerminal[u] = true;
            addEdge(u, parent[u], parentWeight[u]);
            u = parent[u];
        }
    };

    UnionFind regions(z);
    count numberOfComponents = terminals.size();
    for (const Bridge &bridge : bridges) {
        const node regionU = regions.find(closestTerminal[bridge.u]);
        const node regionV = regions.find(closestTerminal[bridge.v]);
        if (regionU != regionV) {
            regions.merge(regionU, regionV);
            --numberOfComponents;
            addEdge(bridge.u, bridge.v, bridge.weight);
            addPathToTerminal(bridge.u);
            addPathToTerminal(bridge.v);
        }
    }

    connected = numberOfComponents <= 1;
    if (!connected) {
        totalWeight = 0;
        tree.removeAllEdges();
    }

    hasRun = true;
}

} /* namespace NetworKit */
